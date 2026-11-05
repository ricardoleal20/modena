//! Bounded framed IO. The calling thread owns and reaps the process; dedicated
//! IO threads keep a blocked pipe write/read inside the request deadline.
use std::io::{Read, Write};
use std::process::{Child, Command, Stdio};
use std::sync::{mpsc, Arc, Mutex};
use std::time::{Duration, Instant};

use serde::de::DeserializeOwned;
use serde::Deserialize;

use super::{failure, Result, MAX_STDERR, MAX_STDOUT, WORKER};

const SESSION_PROTOCOL: &str = "modena-ase-emt-session-v1";
enum Event {
    Written,
    Reply(Vec<u8>),
    Error(String),
}
#[derive(Default)]
struct Stderr {
    bytes: u64,
    diagnostic: Vec<u8>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Frame<T> {
    session_protocol: String,
    request_id: u64,
    response: T,
}

pub(super) struct Session {
    child: Option<Child>,
    input: Option<mpsc::SyncSender<Vec<u8>>>,
    events: mpsc::Receiver<Event>,
    stderr: Arc<Mutex<Stderr>>,
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
        // At most one pending request and two completed IO events. A noisy
        // provider cannot turn frames or diagnostics into an unbounded queue.
        let (input, requests) = mpsc::sync_channel::<Vec<u8>>(1);
        let (events_tx, events) = mpsc::sync_channel(2);
        let writer_events = events_tx.clone();
        std::thread::spawn(move || {
            while let Ok(data) = requests.recv() {
                let result = stdin.write_all(&data).and_then(|()| stdin.flush());
                let event = match result {
                    Ok(()) => Event::Written,
                    Err(e) => Event::Error(format!("session write failed: {e}")),
                };
                let failed = matches!(event, Event::Error(_));
                if writer_events.send(event).is_err() || failed {
                    break;
                }
            }
        });
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
            let event = match reply {
                Ok(data) => Event::Reply(data),
                Err(e) => Event::Error(format!("session read failed (EOF or invalid frame): {e}")),
            };
            let failed = matches!(event, Event::Error(_));
            if events_tx.send(event).is_err() || failed {
                break;
            }
        });
        let stderr = Arc::new(Mutex::new(Stderr::default()));
        let capture = stderr.clone();
        std::thread::spawn(move || {
            let mut chunk = [0; 4096];
            loop {
                match errors.read(&mut chunk) {
                    Ok(0) | Err(_) => break,
                    Ok(n) => {
                        let mut capture = capture.lock().unwrap_or_else(|e| e.into_inner());
                        capture.bytes = capture.bytes.saturating_add(n as u64);
                        let keep = n.min(4096 - capture.diagnostic.len());
                        capture.diagnostic.extend_from_slice(&chunk[..keep]);
                    }
                }
            }
        });
        Ok(Self {
            child: Some(child),
            input: Some(input),
            events,
            stderr,
            next_id: 1,
            timeout: Duration::from_millis(timeout_ms),
        })
    }

    pub(super) fn request<T: DeserializeOwned>(
        &mut self,
        request: &serde_json::Value,
    ) -> Result<T> {
        if self.child.is_none() {
            return Err(failure("persistent session is terminal"));
        }
        let started = Instant::now();
        let result = (|| {
            let id = self.next_id;
            self.next_id = id
                .checked_add(1)
                .ok_or_else(|| failure("session request ID exhausted"))?;
            let data = serde_json::to_vec(&serde_json::json!({"session_protocol":SESSION_PROTOCOL,"request_id":id,"request":request}))
                .map_err(|e| failure(format!("invalid session request: {e}")))?;
            if data.len() as u64 > MAX_STDOUT {
                return Err(failure("session input frame limit exceeded"));
            }
            let mut framed = Vec::with_capacity(data.len() + 4);
            framed.extend_from_slice(&(data.len() as u32).to_be_bytes());
            framed.extend_from_slice(&data);
            self.input
                .as_ref()
                .expect("live session input")
                .try_send(framed)
                .map_err(|e| failure(format!("session input unavailable: {e}")))?;
            let mut written = false;
            let mut reply = None;
            loop {
                if self.stderr.lock().unwrap_or_else(|e| e.into_inner()).bytes > MAX_STDERR {
                    return Err(failure("session stderr output limit exceeded"));
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
                match self.events.recv_timeout(
                    Duration::from_millis(5).min(self.timeout.saturating_sub(started.elapsed())),
                ) {
                    Ok(Event::Written) => written = true,
                    Ok(Event::Reply(data)) if reply.is_none() => reply = Some(data),
                    Ok(Event::Reply(_)) => return Err(failure("unsolicited session reply")),
                    Ok(Event::Error(e)) => return Err(failure(e)),
                    Err(mpsc::RecvTimeoutError::Timeout) => (),
                    Err(e) => return Err(failure(format!("session IO disconnected: {e}"))),
                }
            }
            let frame: Frame<T> = serde_json::from_slice(&reply.expect("received reply"))
                .map_err(|e| failure(format!("invalid session response: {e}")))?;
            if frame.session_protocol != SESSION_PROTOCOL || frame.request_id != id {
                return Err(failure("session protocol or request ID mismatch"));
            }
            Ok(frame.response)
        })();
        if let Err(e) = result {
            self.terminate();
            let capture = self.stderr.lock().unwrap_or_else(|e| e.into_inner());
            return Err(failure(format!(
                "{e}; stderr: {}",
                String::from_utf8_lossy(&capture.diagnostic)
            )));
        }
        result
    }

    pub(super) fn terminate(&mut self) {
        if let Some(mut child) = self.child.take() {
            let _ = child.kill();
            let _ = child.wait();
        }
        self.input.take();
    }
}
impl Drop for Session {
    fn drop(&mut self) {
        self.terminate();
    }
}
