//! Bounded framed IO with independently enforced process lifetime.
use std::io::{Read, Write};
use std::process::{Child, Command, Stdio};
use std::sync::{mpsc, Arc, Mutex, Weak};
use std::time::{Duration, Instant};

use rand::TryRng;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};

use super::{failure, Result, MAX_STDERR, MAX_STDOUT, WORKER};

const SESSION_PROTOCOL: &str = "modena-ase-emt-session-v1";
enum Event {
    Written,
    Reply(Vec<u8>),
}
#[derive(Default)]
struct Stderr {
    bytes: u64,
    diagnostic: Vec<u8>,
}
struct Process {
    child: Option<Child>,
    failure: Option<String>,
}
struct Shared {
    process: Mutex<Process>,
    input: Mutex<Option<mpsc::SyncSender<Vec<u8>>>>,
    events: Mutex<Option<mpsc::Receiver<Event>>>,
    pending: Mutex<Option<(u64, [u8; 32])>>,
    stderr: Mutex<Stderr>,
}
impl Shared {
    fn terminate(&self, reason: &str) {
        // Keep state locked through wait: observing terminal state also
        // observes completed reaping, even when another IO thread failed.
        {
            let mut process = self.process.lock().unwrap_or_else(|e| e.into_inner());
            process
                .failure
                .get_or_insert_with(|| failure(reason).to_string());
            if let Some(mut child) = process.child.take() {
                let _ = child.kill();
                let _ = child.wait();
            }
        }
        self.input.lock().unwrap_or_else(|e| e.into_inner()).take();
        // Disconnect producers and release queued frames even if the failed
        // provider object remains alive indefinitely.
        self.events.lock().unwrap_or_else(|e| e.into_inner()).take();
        self.pending
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .take();
    }
    fn failure(&self) -> Option<String> {
        self.process
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .failure
            .clone()
    }
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Frame<T> {
    session_protocol: String,
    request_id: u64,
    nonce: [u8; 32],
    response: T,
}

// Header and JSON share one allocation, without intermediate Value trees or
// a second framed copy. The cap is applied during every serializer write.
struct CappedFrame(Vec<u8>);
impl Write for CappedFrame {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let limit = MAX_STDOUT as usize + 4;
        if bytes.len() > limit - self.0.len() {
            return Err(std::io::Error::other("session input frame limit exceeded"));
        }
        let required = self.0.len() + bytes.len();
        if required > self.0.capacity() {
            let capacity = required.max(self.0.capacity().saturating_mul(2)).min(limit);
            self.0.reserve_exact(capacity - self.0.len());
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encode_frame(id: u64, nonce: &[u8; 32], request: &impl Serialize) -> Result<Vec<u8>> {
    #[derive(Serialize)]
    struct RequestFrame<'a, T> {
        session_protocol: &'static str,
        request_id: u64,
        nonce: &'a [u8; 32],
        request: &'a T,
    }
    let mut framed = CappedFrame(Vec::with_capacity(4096));
    framed.0.extend_from_slice(&[0; 4]);
    serde_json::to_writer(
        &mut framed,
        &RequestFrame {
            session_protocol: SESSION_PROTOCOL,
            request_id: id,
            nonce,
            request,
        },
    )
    .map_err(|e| failure(format!("cannot serialize session request: {e}")))?;
    let size = (framed.0.len() - 4) as u32;
    framed.0[..4].copy_from_slice(&size.to_be_bytes());
    Ok(framed.0)
}
fn send_event(shared: &Weak<Shared>, tx: &mpsc::SyncSender<Event>, event: Event) -> bool {
    // No producer remains blocked behind an abandoned/full event queue.
    if tx.try_send(event).is_err() {
        if let Some(shared) = shared.upgrade() {
            shared.terminate("unsolicited session output or disconnected IO queue");
        }
        false
    } else {
        true
    }
}

pub(super) struct Session {
    shared: Arc<Shared>,
    next_id: u64,
    timeout: Duration,
}
impl Session {
    pub(super) fn new(python: &str, timeout_ms: u64) -> Result<Self> {
        if python.trim().is_empty() || timeout_ms == 0 {
            return Err(failure("python and positive timeout_ms required"));
        }
        let mut child = Command::new(python)
            .args(["-u", "-c", WORKER, "--persistent"])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| failure(format!("cannot launch Python: {e}")))?;
        let mut stdin = child.stdin.take().expect("piped stdin");
        let mut stdout = child.stdout.take().expect("piped stdout");
        let mut errors = child.stderr.take().expect("piped stderr");
        let (input, requests) = mpsc::sync_channel::<Vec<u8>>(1);
        let (events_tx, events) = mpsc::sync_channel(2);
        let shared = Arc::new(Shared {
            process: Mutex::new(Process {
                child: Some(child),
                failure: None,
            }),
            input: Mutex::new(Some(input)),
            events: Mutex::new(Some(events)),
            pending: Mutex::new(None),
            stderr: Mutex::new(Stderr::default()),
        });
        // Threads keep only Weak references while waiting on pipes. IO
        // failures can reap idle workers without extending provider lifetime.
        let writer_shared = Arc::downgrade(&shared);
        let writer_events = events_tx.clone();
        std::thread::spawn(move || {
            while let Ok(data) = requests.recv() {
                if let Err(e) = stdin.write_all(&data).and_then(|()| stdin.flush()) {
                    if let Some(shared) = writer_shared.upgrade() {
                        shared.terminate(&format!("session write failed: {e}"));
                    }
                    break;
                }
                if !send_event(&writer_shared, &writer_events, Event::Written) {
                    break;
                }
            }
        });
        let reader_shared = Arc::downgrade(&shared);
        std::thread::spawn(move || loop {
            let reply = (|| -> std::io::Result<Vec<u8>> {
                let mut header = [0; 4];
                stdout.read_exact(&mut header)?;
                let size = u32::from_be_bytes(header) as u64;
                if size == 0 || size > MAX_STDOUT {
                    return Err(std::io::Error::other("session output frame limit exceeded"));
                }
                let mut data = vec![0; size as usize];
                stdout.read_exact(&mut data)?;
                Ok(data)
            })();
            match reply {
                Ok(data) => {
                    // Validate correlation immediately, including while the
                    // provider is idle. Ignore property contents here rather
                    // than allocating an intermediate response Value tree.
                    let checked = (|| -> Result<()> {
                        let envelope: Frame<serde::de::IgnoredAny> = serde_json::from_slice(&data)
                            .map_err(|e| failure(format!("invalid session envelope: {e}")))?;
                        let Some(shared) = reader_shared.upgrade() else {
                            return Err(failure("session dropped"));
                        };
                        let mut pending = shared.pending.lock().unwrap_or_else(|e| e.into_inner());
                        if envelope.session_protocol != SESSION_PROTOCOL
                            || pending.as_ref() != Some(&(envelope.request_id, envelope.nonce))
                        {
                            return Err(failure(
                                "unsolicited session reply or correlation mismatch",
                            ));
                        }
                        pending.take();
                        Ok(())
                    })();
                    if let Err(e) = checked {
                        if let Some(shared) = reader_shared.upgrade() {
                            shared.terminate(&e.to_string());
                        }
                        break;
                    }
                    if !send_event(&reader_shared, &events_tx, Event::Reply(data)) {
                        break;
                    }
                }
                Err(e) => {
                    if let Some(shared) = reader_shared.upgrade() {
                        shared
                            .terminate(&format!("session read failed (EOF or invalid frame): {e}"));
                    }
                    break;
                }
            }
        });
        let stderr_shared = Arc::downgrade(&shared);
        std::thread::spawn(move || {
            let mut chunk = [0; 4096];
            loop {
                match errors.read(&mut chunk) {
                    Ok(0) => break,
                    Err(e) => {
                        if let Some(shared) = stderr_shared.upgrade() {
                            shared.terminate(&format!("session stderr read failed: {e}"));
                        }
                        break;
                    }
                    Ok(n) => {
                        let Some(shared) = stderr_shared.upgrade() else {
                            break;
                        };
                        let overflow = {
                            let mut capture =
                                shared.stderr.lock().unwrap_or_else(|e| e.into_inner());
                            capture.bytes = capture.bytes.saturating_add(n as u64);
                            let keep = n.min(4096 - capture.diagnostic.len());
                            capture.diagnostic.extend_from_slice(&chunk[..keep]);
                            capture.bytes > MAX_STDERR
                        };
                        if overflow {
                            shared.terminate("session stderr output limit exceeded");
                            break;
                        }
                    }
                }
            }
        });
        Ok(Self {
            shared,
            next_id: 1,
            timeout: Duration::from_millis(timeout_ms),
        })
    }
    pub(super) fn request<T: DeserializeOwned>(&mut self, request: &impl Serialize) -> Result<T> {
        if let Some(reason) = self.shared.failure() {
            return Err(failure(format!("persistent session is terminal: {reason}")));
        }
        let started = Instant::now();
        let result = (|| {
            let id = self.next_id;
            self.next_id = id
                .checked_add(1)
                .ok_or_else(|| failure("session request ID exhausted"))?;
            // Fresh OS entropy prevents buffered replies for predictable
            // integer IDs from impersonating responses to future snapshots.
            let mut nonce = [0; 32];
            rand::rngs::SysRng
                .try_fill_bytes(&mut nonce)
                .map_err(|e| failure(format!("cannot create session nonce: {e}")))?;
            let framed = encode_frame(id, &nonce, request)?;
            *self
                .shared
                .pending
                .lock()
                .unwrap_or_else(|e| e.into_inner()) = Some((id, nonce));
            {
                let input = self.shared.input.lock().unwrap_or_else(|e| e.into_inner());
                input
                    .as_ref()
                    .ok_or_else(|| failure("persistent session is terminal"))?
                    .try_send(framed)
                    .map_err(|e| failure(format!("session input unavailable: {e}")))?;
            }
            let mut written = false;
            let mut reply = None;
            loop {
                if let Some(reason) = self.shared.failure() {
                    return Err(failure(format!("persistent session is terminal: {reason}")));
                }
                if started.elapsed() >= self.timeout {
                    return Err(failure(format!(
                        "session timeout after {} ms",
                        self.timeout.as_millis()
                    )));
                }
                if written && reply.is_some() {
                    break;
                }
                let event = {
                    let events = self.shared.events.lock().unwrap_or_else(|e| e.into_inner());
                    events
                        .as_ref()
                        .ok_or_else(|| failure("persistent session is terminal"))?
                        .recv_timeout(
                            Duration::from_millis(5)
                                .min(self.timeout.saturating_sub(started.elapsed())),
                        )
                };
                match event {
                    Ok(Event::Written) => written = true,
                    Ok(Event::Reply(data)) if reply.is_none() => reply = Some(data),
                    Ok(Event::Reply(_)) => return Err(failure("unsolicited session reply")),
                    Err(mpsc::RecvTimeoutError::Timeout) => (),
                    Err(e) => return Err(failure(format!("session IO disconnected: {e}"))),
                }
            }
            let frame: Frame<T> = serde_json::from_slice(&reply.expect("received reply"))
                .map_err(|e| failure(format!("invalid session response: {e}")))?;
            if frame.session_protocol != SESSION_PROTOCOL
                || frame.request_id != id
                || frame.nonce != nonce
            {
                return Err(failure("session protocol, request ID or nonce mismatch"));
            }
            Ok(frame.response)
        })();
        if let Err(e) = result {
            self.shared
                .terminate("request or response validation failed");
            let capture = self.shared.stderr.lock().unwrap_or_else(|e| e.into_inner());
            return Err(failure(format!(
                "{e}; stderr: {}",
                String::from_utf8_lossy(&capture.diagnostic)
            )));
        }
        result
    }
    pub(super) fn terminate(&mut self) {
        self.shared
            .terminate("provider validation failed or dropped");
    }
}
impl Drop for Session {
    fn drop(&mut self) {
        self.terminate();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn terminalization_disconnects_a_full_io_queue() {
        let (tx, events) = mpsc::sync_channel(1);
        assert!(tx.send(Event::Reply(vec![0; 1024])).is_ok());
        let (input, _requests) = mpsc::sync_channel(1);
        let shared = Shared {
            process: Mutex::new(Process {
                child: None,
                failure: None,
            }),
            input: Mutex::new(Some(input)),
            events: Mutex::new(Some(events)),
            pending: Mutex::new(None),
            stderr: Mutex::new(Stderr::default()),
        };
        shared.terminate("test failure");
        assert!(
            matches!(
                tx.try_send(Event::Written),
                Err(mpsc::TrySendError::Disconnected(_))
            ),
            "terminal session retained its full IO queue"
        );
        assert!(shared.events.lock().unwrap().is_none());
    }
    #[test]
    fn oversized_serialization_stops_before_materializing_the_entire_request() {
        use serde::ser::SerializeSeq;
        use std::sync::atomic::{AtomicUsize, Ordering};
        struct Growing(AtomicUsize);
        impl Serialize for Growing {
            fn serialize<S: serde::Serializer>(
                &self,
                serializer: S,
            ) -> std::result::Result<S::Ok, S::Error> {
                let chunk = "x".repeat(4096);
                let mut seq = serializer.serialize_seq(Some(8192))?;
                for _ in 0..8192 {
                    self.0.fetch_add(1, Ordering::Relaxed);
                    seq.serialize_element(&chunk)?;
                }
                seq.end()
            }
        }
        let payload = Growing(AtomicUsize::new(0));
        assert!(encode_frame(1, &[0; 32], &payload).is_err());
        assert!(
            payload.0.load(Ordering::Relaxed) < 8192,
            "serialized the complete oversized request before enforcing the cap"
        );
        let mut writer = CappedFrame(Vec::with_capacity(4096));
        writer.0.resize(4, 0);
        let chunk = vec![0; 4096];
        while writer.write_all(&chunk).is_ok() {}
        assert!(writer.0.capacity() <= MAX_STDOUT as usize + 4);
        assert!(writer.0.len() <= MAX_STDOUT as usize + 4);
    }
}
