//! Conservative process-wide requested Rust heap cap, including opaque core/FFT storage.
//! This is neither a per-model capacity count nor allocator overhead/OS RSS.
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

pub const PROCESS_HEAP_CAP: usize = 268_435_456;
pub const HEAP_EXIT: i32 = 86;
pub const TIMING_EXIT: i32 = 87;
pub const SYSTEM_ALLOC_EXIT: i32 = 88;
pub const SYSTEM_ALLOC_MARKER: &[u8] = b"MODENA_B05_SYSTEM_ALLOCATION_FAILED\n";
pub const HEAP_MARKER: &[u8] = b"MODENA_B05_RUST_HEAP_CAP\n";
pub const TIMING_MARKER: &[u8] = b"MODENA_B05_TIMING_ALLOCATOR_CALLBACK\n";

/// Reservations include in-flight allocations. Failed system allocations roll back.
/// SeqCst updates make the process cap conservative across Rayon workers.
pub struct HeapCounter {
    live: AtomicUsize,
    peak: AtomicUsize,
    cap: usize,
}
impl HeapCounter {
    pub const fn new(cap: usize) -> Self {
        Self {
            live: AtomicUsize::new(0),
            peak: AtomicUsize::new(0),
            cap,
        }
    }
    pub fn reserve(&self, bytes: usize) -> bool {
        let mut old = self.live.load(Ordering::SeqCst);
        loop {
            let Some(next) = old.checked_add(bytes).filter(|n| *n <= self.cap) else {
                return false;
            };
            match self
                .live
                .compare_exchange_weak(old, next, Ordering::SeqCst, Ordering::SeqCst)
            {
                Ok(_) => {
                    self.peak.fetch_max(next, Ordering::SeqCst);
                    return true;
                }
                Err(current) => old = current,
            }
        }
    }
    pub fn release(&self, bytes: usize) -> bool {
        self.live
            .fetch_update(Ordering::SeqCst, Ordering::SeqCst, |live| {
                live.checked_sub(bytes)
            })
            .is_ok()
    }
    pub fn live(&self) -> usize {
        self.live.load(Ordering::SeqCst)
    }
    pub fn peak(&self) -> usize {
        self.peak.load(Ordering::SeqCst)
    }
}
static HEAP: HeapCounter = HeapCounter::new(PROCESS_HEAP_CAP);
/// SeqCst admission state shared by callbacks and guard entry.
pub struct CallbackState {
    active: AtomicUsize,
    forbidden: AtomicBool,
}
impl CallbackState {
    pub const fn new() -> Self {
        Self {
            active: AtomicUsize::new(0),
            forbidden: AtomicBool::new(false),
        }
    }
    pub fn enter_callback(&self) -> bool {
        self.active
            .fetch_update(Ordering::SeqCst, Ordering::SeqCst, |n| n.checked_add(1))
            .is_ok()
            && !self.forbidden.load(Ordering::SeqCst)
    }
    pub fn leave_callback(&self) {
        let _ = self
            .active
            .fetch_update(Ordering::SeqCst, Ordering::SeqCst, |n| n.checked_sub(1));
    }
    pub fn begin_guard(&self) -> bool {
        if self
            .forbidden
            .compare_exchange(false, true, Ordering::SeqCst, Ordering::SeqCst)
            .is_err()
        {
            return false;
        }
        // SeqCst proof: a callback increments active BEFORE reading forbidden.
        // If it read the old false, its increment precedes our CAS and is seen
        // here unless its entire allocator body already finished. A later
        // increment must read true and fail closed. Thus no accepted guard can
        // miss a callback that observed the old flag. Both admission and this
        // quiescence check occur before the caller starts any clock.
        if self.active.load(Ordering::SeqCst) != 0 {
            self.forbidden.store(false, Ordering::SeqCst);
            return false;
        }
        true
    }
    pub fn end_guard(&self) {
        self.forbidden.store(false, Ordering::SeqCst);
    }
}
impl Default for CallbackState {
    fn default() -> Self {
        Self::new()
    }
}
static CALLBACK_STATE: CallbackState = CallbackState::new();

/// Enter/exit OUTSIDE primary timing clocks. Any allocator callback fails closed,
/// including deallocation and realloc shrink. Use in an isolated campaign process.
#[must_use]
pub struct NoAllocationGuard {
    _private: (),
}
impl NoAllocationGuard {
    pub fn enter() -> Option<Self> {
        CALLBACK_STATE.begin_guard().then(|| Self { _private: () })
    }
}
impl Drop for NoAllocationGuard {
    fn drop(&mut self) {
        CALLBACK_STATE.end_guard();
    }
}
pub fn requested_heap() -> usize {
    HEAP.live()
}
pub fn peak_requested_heap() -> usize {
    HEAP.peak()
}

// The supported first campaign host is Unix/macOS. Avoid stdio, formatting,
// unwinding, destructors and exit hooks here: each could reenter the allocator.
#[cfg(unix)]
fn terminate(marker: &[u8], code: i32) -> ! {
    unsafe extern "C" {
        fn write(fd: i32, buf: *const u8, count: usize) -> isize;
        fn _exit(status: i32) -> !;
    }
    let mut sent = 0;
    while sent < marker.len() {
        // SAFETY: marker is a static valid byte slice; no buffer survives write.
        let n = unsafe { write(2, marker.as_ptr().add(sent), marker.len() - sent) };
        if n <= 0 {
            break;
        }
        sent += n as usize;
    }
    // SAFETY: immediate process termination is the deliberate hard-cap contract.
    unsafe { _exit(code) }
}
#[cfg(not(unix))]
compile_error!("B05 capped allocator requires the registered Unix campaign host");
fn callback(marker: &[u8]) {
    if !CALLBACK_STATE.enter_callback() {
        terminate(marker, TIMING_EXIT);
    }
}
fn reserve(bytes: usize) {
    if !HEAP.reserve(bytes) {
        terminate(HEAP_MARKER, HEAP_EXIT);
    }
}
pub struct CappedAllocator;
// SAFETY: all pointer/layout operations delegate unchanged to System. Atomic
// reservations are taken BEFORE allocation/growth, rolled back on null, and
// released only after System deallocation/shrink. No callback unwinds or allocates.
unsafe impl GlobalAlloc for CappedAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        callback(b"MODENA_B05_TIMING_ALLOCATOR_CALLBACK\nkind=alloc\n");
        reserve(layout.size());
        // SAFETY: caller supplies the valid layout required by GlobalAlloc.
        let ptr = unsafe { System.alloc(layout) };
        if ptr.is_null() {
            HEAP.release(layout.size());
            terminate(SYSTEM_ALLOC_MARKER, SYSTEM_ALLOC_EXIT);
        }
        CALLBACK_STATE.leave_callback();
        ptr
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        callback(b"MODENA_B05_TIMING_ALLOCATOR_CALLBACK\nkind=alloc_zeroed\n");
        reserve(layout.size());
        // SAFETY: caller supplies the valid layout required by GlobalAlloc.
        let ptr = unsafe { System.alloc_zeroed(layout) };
        if ptr.is_null() {
            HEAP.release(layout.size());
            terminate(SYSTEM_ALLOC_MARKER, SYSTEM_ALLOC_EXIT);
        }
        CALLBACK_STATE.leave_callback();
        ptr
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        callback(b"MODENA_B05_TIMING_ALLOCATOR_CALLBACK\nkind=dealloc\n");
        // SAFETY: caller supplies the live System allocation and original layout.
        unsafe {
            System.dealloc(ptr, layout);
        }
        HEAP.release(layout.size());
        CALLBACK_STATE.leave_callback();
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        callback(b"MODENA_B05_TIMING_ALLOCATOR_CALLBACK\nkind=realloc\n");
        let growth = new_size.saturating_sub(layout.size());
        reserve(growth);
        // SAFETY: GlobalAlloc requires valid ptr/layout and nonzero new_size.
        let next = unsafe { System.realloc(ptr, layout, new_size) };
        if next.is_null() {
            HEAP.release(growth);
            terminate(SYSTEM_ALLOC_MARKER, SYSTEM_ALLOC_EXIT);
        } else if new_size < layout.size() {
            HEAP.release(layout.size() - new_size);
        }
        CALLBACK_STATE.leave_callback();
        next
    }
}
#[global_allocator]
static ALLOCATOR: CappedAllocator = CappedAllocator;
