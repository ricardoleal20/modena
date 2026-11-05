use std::alloc::{GlobalAlloc, Layout, System as AllocSystem};
use std::sync::atomic::{AtomicUsize, Ordering};

pub struct CountingAllocator;
pub static FORCE_BYTES: AtomicUsize = AtomicUsize::new(0);
pub static FORCE_ALLOCS: AtomicUsize = AtomicUsize::new(0);
pub static FORCE_REALLOCS: AtomicUsize = AtomicUsize::new(0);
pub static ALLOCATED_BYTES: AtomicUsize = AtomicUsize::new(0);
pub static LIVE_BYTES: AtomicUsize = AtomicUsize::new(0);
pub static HIGH_WATER: AtomicUsize = AtomicUsize::new(0);

fn allocated(bytes: usize, realloc: bool) {
    let live = LIVE_BYTES.fetch_add(bytes, Ordering::SeqCst) + bytes;
    HIGH_WATER.fetch_max(live, Ordering::SeqCst);
    if bytes == FORCE_BYTES.load(Ordering::SeqCst) {
        if realloc {
            FORCE_REALLOCS.fetch_add(1, Ordering::SeqCst);
        } else {
            FORCE_ALLOCS.fetch_add(1, Ordering::SeqCst);
        }
        ALLOCATED_BYTES.fetch_add(bytes, Ordering::SeqCst);
    }
}
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { AllocSystem.alloc(layout) };
        if !ptr.is_null() {
            allocated(layout.size(), false);
        }
        ptr
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { AllocSystem.alloc_zeroed(layout) };
        if !ptr.is_null() {
            allocated(layout.size(), false);
        }
        ptr
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        LIVE_BYTES.fetch_sub(layout.size(), Ordering::SeqCst);
        unsafe { AllocSystem.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        let new = unsafe { AllocSystem.realloc(ptr, layout, size) };
        if !new.is_null() {
            LIVE_BYTES.fetch_sub(layout.size(), Ordering::SeqCst);
            allocated(size, true);
        }
        new
    }
}

pub fn reset(bytes: usize) {
    FORCE_BYTES.store(bytes, Ordering::SeqCst);
    FORCE_ALLOCS.store(0, Ordering::SeqCst);
    FORCE_REALLOCS.store(0, Ordering::SeqCst);
    ALLOCATED_BYTES.store(0, Ordering::SeqCst);
    HIGH_WATER.store(LIVE_BYTES.load(Ordering::SeqCst), Ordering::SeqCst);
}

pub fn replicas(n: usize) -> [u32; 3] {
    match n {
        864 => [6, 6, 6],
        2048 => [8, 8, 8],
        5000 => [10, 5, 25],
        10000 => [10, 10, 25],
        _ => panic!("supported sizes: 864, 2048, 5000, 10000"),
    }
}

pub fn system(n: usize) -> modena_core::System {
    modena_core::lattice::fcc(1.6796, replicas(n), "Ar")
}
