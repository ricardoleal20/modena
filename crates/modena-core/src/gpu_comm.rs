//! GPU offload of the Lennard-Jones pair kernel (M9, `gpu` feature).
//!
//! A [`ForceModel`] backed by WGSL compute shaders running on the GPU via
//! `wgpu`, evaluating the *same* shifted-force Lennard-Jones potential as the
//! CPU path, with orthorhombic minimum-image periodicity. Two gather
//! strategies ([`GpuGather`]):
//!
//! - [`GpuGather::AllPairs`] — the original O(N²) kernel: one thread per
//!   atom `i`, looping over *all* `j`. Every pair is visited twice (once from
//!   each end), so the per-atom energy/virial are halved and summed on the
//!   host. Kept as the reference/**determinism oracle** for the cell-list
//!   kernel (tests compare the two, and against the serial f64 path).
//! - [`GpuGather::CellList`] — the G3 remedy for the all-pairs kernel's
//!   gather bottleneck (measured: the j-loop gather + minimum-image + distance
//!   work is ~85–90% of the all-pairs kernel time, because only ~55 of N
//!   candidates per atom are within the cutoff). The CPU cell-list/Verlet
//!   bins ([`NeighborList`]) are converted into contiguous per-atom j-ranges
//!   (a CSR: `j_start[i]..j_start[i+1]` into a flat `j_list`) and uploaded;
//!   the GPU then runs one thread per `i` over just its own range — no
//!   atomics, no per-thread O(N) scan. Each pair `(i, j)` appears once in
//!   `i`'s range and once in `j`'s, so per-atom energy/virial are halved
//!   (the host sum then counts each pair exactly once) while each atom's
//!   force is its own full net force.
//!
//! ## Precision
//!
//! The shaders compute in `f32`; the CPU path is `f64`. The two therefore agree
//! to single precision (~1e-5 relative), not bit-for-bit. This is the standard
//! GPU MD trade-off (LAMMPS's GPU package does the same). Accumulation uses
//! Kahan compensation in-f32, so the summation error is O(1) instead of
//! O(n·eps) over the pair loop; the host reduces f32 outputs into f64
//! (mixed-mode style). A double-precision shader is **not possible on Apple
//! silicon** (WGSL has no f64 over Metal — gpuweb#2805), so the f32 tier is the
//! shipped form, validated as a separate single-precision path in
//! `tests/gpu_pair.rs`.
//!
//! ## Scope
//!
//! Single-species (reduced) LJ on an orthorhombic box — the Phase-1 benchmark
//! potential. Uploading a per-type-pair `(σ, ε)` table + a triclinic minimum-
//! image path are straightforward extensions of the same buffer layout. The
//! seam is the [`ForceModel`] trait: the runner treats the GPU model exactly
//! like any other.
//!
//! ## Graceful absence
//!
//! [`GpuPairForceModel::new`] returns `None` when no suitable GPU adapter is
//! available (headless CI, no drivers) or a shader fails to compile, so
//! callers and tests skip cleanly instead of panicking.

#![cfg(feature = "gpu")]

use crate::error::{ModenaError, Result};
use crate::force_model::{ForceModel, Virial};
use crate::neighbor::NeighborList;
use crate::system::System;

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;
use wgpu::PollType;

/// Uniform parameters, packed as three `vec4<f32>` (48 bytes) to match the
/// WGSL `array<vec4<f32>, 3>` layout with no alignment surprises.
///
/// `[n, sigma, epsilon, cutoff]`, `[cutoff2, u_rc, f_rc, _]`, `[bx, by, bz, _]`.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct GpuParams {
    n: f32,
    sigma: f32,
    epsilon: f32,
    cutoff: f32,
    cutoff2: f32,
    u_rc: f32,
    f_rc: f32,
    _pad0: f32,
    bx: f32,
    by: f32,
    bz: f32,
    _pad1: f32,
}

const SHADER: &str = r#"
@group(0) @binding(0) var<uniform> params: array<vec4<f32>, 3>;
@group(0) @binding(1) var<storage, read>       positions: array<f32>;
@group(0) @binding(2) var<storage, read_write> forces:   array<f32>;
@group(0) @binding(3) var<storage, read_write> aux:      array<f32>; // [pe_i (n), virial_i (n)]

fn min_image(d: vec3<f32>, box: vec3<f32>) -> vec3<f32> {
    return d - box * round(d / box);
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n      = u32(params[0].x);
    let sigma  = params[0].y;
    let eps    = params[0].z;
    let cutoff = params[0].w;
    let cutoff2 = params[1].x;
    let u_rc   = params[1].y;
    let f_rc   = params[1].z;
    let box    = vec3<f32>(params[2].x, params[2].y, params[2].z);

    let i = gid.x;
    if (i >= n) {
        return;
    }
    let ri = vec3<f32>(positions[3u*i], positions[3u*i + 1u], positions[3u*i + 2u]);

    var f = vec3<f32>(0.0);
    var f_comp = vec3<f32>(0.0);
    var pe = 0.0;
    var pe_comp = 0.0;
    var vir = 0.0;
    var vir_comp = 0.0;

    let sigma2 = sigma * sigma;
    for (var j = 0u; j < n; j = j + 1u) {
        if (j == i) {
            continue;
        }
        let rj = vec3<f32>(positions[3u*j], positions[3u*j + 1u], positions[3u*j + 2u]);
        let d = min_image(ri - rj, box);
        let r2 = dot(d, d);
        if ((r2 >= cutoff2) || (r2 < 1e-14)) {
            continue;
        }
        let r = sqrt(r2);
        let sr2 = sigma2 / r2;
        let sr6 = sr2 * sr2 * sr2;
        let sr12 = sr6 * sr6;
        let u_raw = 4.0 * eps * (sr12 - sr6);
        let f_raw = (24.0 * eps / r) * (2.0 * sr12 - sr6);
        // Shifted-force form, matching the CPU `LennardJones` kernel.
        let u = u_raw - u_rc + (r - cutoff) * f_rc;
        let f_over_r = (f_raw - f_rc) / r;
        let fv = d * f_over_r;
        let u2 = 0.5 * u;        // each pair counted from both ends
        let v2 = 0.5 * dot(d, fv);
        // Kahan (compensated) accumulation: f32 registers, O(1) error instead
        // of O(n·eps) over the all-pairs loop. Host-side f64 reduction still
        // applies (LAMMPS mixed-mode style).
        let fy = fv - f_comp;
        let ft = f + fy;
        f_comp = (ft - f) - fy;
        f = ft;
        let pey = u2 - pe_comp;
        let pet = pe + pey;
        pe_comp = (pet - pe) - pey;
        pe = pet;
        let viry = v2 - vir_comp;
        let virt = vir + viry;
        vir_comp = (virt - vir) - viry;
        vir = virt;
    }

    forces[3u*i]      = f.x;
    forces[3u*i + 1u] = f.y;
    forces[3u*i + 2u] = f.z;
    aux[i]     = pe;
    aux[n + i] = vir;
}
"#;

/// Cell-list (CSR) kernel: one thread per atom `i`, iterating only the
/// per-atom j-range `[j_start[i], j_start[i+1])` built on the CPU from the
/// Verlet-binned half list. Each pair `(i, j)` appears once in `i`'s range and
/// once in `j`'s, so per-atom energy/virial are halved (the host sum then
/// counts each pair exactly once) while each atom's force is its own full net
/// force. Same shifted-force LJ + Kahan f32 accumulation as the all-pairs
/// kernel.
const SHADER_CELL_LIST: &str = r#"
@group(0) @binding(0) var<uniform> params: array<vec4<f32>, 3>;
@group(0) @binding(1) var<storage, read>       positions: array<f32>;
@group(0) @binding(2) var<storage, read_write> forces:   array<f32>;
@group(0) @binding(3) var<storage, read_write> aux:      array<f32>;
@group(0) @binding(4) var<storage, read>       j_start:  array<u32>; // n + 1 entries
@group(0) @binding(5) var<storage, read>       j_list:   array<u32>;

fn min_image(d: vec3<f32>, box: vec3<f32>) -> vec3<f32> {
    return d - box * round(d / box);
}

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n      = u32(params[0].x);
    let sigma  = params[0].y;
    let eps    = params[0].z;
    let cutoff = params[0].w;
    let cutoff2 = params[1].x;
    let u_rc   = params[1].y;
    let f_rc   = params[1].z;
    let box    = vec3<f32>(params[2].x, params[2].y, params[2].z);

    let i = gid.x;
    if (i >= n) {
        return;
    }
    let ri = vec3<f32>(positions[3u*i], positions[3u*i + 1u], positions[3u*i + 2u]);

    var f = vec3<f32>(0.0);
    var f_comp = vec3<f32>(0.0);
    var pe = 0.0;
    var pe_comp = 0.0;
    var vir = 0.0;
    var vir_comp = 0.0;

    let sigma2 = sigma * sigma;
    let begin = j_start[i];
    let end   = j_start[i + 1u];
    for (var k = begin; k < end; k = k + 1u) {
        let j = j_list[k];
        let rj = vec3<f32>(positions[3u*j], positions[3u*j + 1u], positions[3u*j + 2u]);
        let d = min_image(ri - rj, box);
        let r2 = dot(d, d);
        if ((r2 >= cutoff2) || (r2 < 1e-14)) {
            continue;
        }
        let r = sqrt(r2);
        let sr2 = sigma2 / r2;
        let sr6 = sr2 * sr2 * sr2;
        let sr12 = sr6 * sr6;
        let u_raw = 4.0 * eps * (sr12 - sr6);
        let f_raw = (24.0 * eps / r) * (2.0 * sr12 - sr6);
        // Shifted-force form, matching the CPU `LennardJones` kernel.
        let u = u_raw - u_rc + (r - cutoff) * f_rc;
        let f_over_r = (f_raw - f_rc) / r;
        let fv = d * f_over_r;
        // Halved per-atom energy/virial: every pair appears in *both*
        // endpoints' j-ranges (the half list feeds i's and j's ranges), so the
        // host sum would count each pair twice without the 0.5. Forces are
        // NOT halved — each atom needs its own full net force.
        let u2 = 0.5 * u;
        let v2 = 0.5 * dot(d, fv);
        // Kahan (compensated) accumulation: f32 registers, O(1) error instead
        // of O(n·eps). Host-side f64 reduction still applies.
        let fy = fv - f_comp;
        let ft = f + fy;
        f_comp = (ft - f) - fy;
        f = ft;
        let pey = u2 - pe_comp;
        let pet = pe + pey;
        pe_comp = (pet - pe) - pey;
        pe = pet;
        let viry = v2 - vir_comp;
        let virt = vir + viry;
        vir_comp = (virt - vir) - viry;
        vir = virt;
    }

    forces[3u*i]      = f.x;
    forces[3u*i + 1u] = f.y;
    forces[3u*i + 2u] = f.z;
    aux[i]     = pe;
    aux[n + i] = vir;
}
"#;

/// Which j-candidate source the GPU pair kernel gathers from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GpuGather {
    /// All-pairs O(N²) j-loop — the original kernel. Every pair is visited
    /// from both ends (per-atom energy/virial halved). Kept as the
    /// reference/determinism oracle for [`GpuGather::CellList`].
    #[default]
    AllPairs,
    /// Cell-list gather (the G3 remedy): CPU Verlet bins → contiguous
    /// per-atom j-ranges (CSR) uploaded to the GPU; one thread per `i`,
    /// no atomics. Requires the caller's [`NeighborList`] (the standard
    /// runner seam already provides it).
    CellList,
}

/// GPU-accelerated Lennard-Jones [`ForceModel`].
///
/// Construct with [`new`](Self::new) (all-pairs reference kernel) or
/// [`with_gather`](Self::with_gather) (e.g. [`GpuGather::CellList`]); both
/// return `None` when no usable GPU adapter exists (so CI without a GPU skips
/// gracefully). The model owns the `wgpu` device/queue and pre-compiled
/// compute pipelines; work buffers + bind groups are cached per atom count
/// (see [`ensure_buffers`](Self::ensure_buffers)) — a previous per-call
/// `create_buffer`×4 cost ~2 ms/step and made the GPU slower than the CPU
/// below ~1k atoms.
pub struct GpuPairForceModel {
    device: wgpu::Device,
    queue: wgpu::Queue,
    /// All-pairs reference pipeline (bindings 0–3).
    pipeline: wgpu::ComputePipeline,
    /// Cell-list (CSR) pipeline (bindings 0–5). Only used with
    /// [`GpuGather::CellList`].
    pipeline_cell_list: wgpu::ComputePipeline,
    /// Layout for the CSR pipeline's extra read-only storage bindings.
    cell_list_bgl: wgpu::BindGroupLayout,
    uniform: wgpu::Buffer,
    params: GpuParams,
    gather: GpuGather,
    /// Buffers + bind groups cached per atom count: creating four buffers and a
    /// bind group per `compute()` call cost ~2 ms/step (the reason the GPU is
    /// slower than the CPU below ~1k atoms). Reallocated only when `n`
    /// changes.
    buffers: std::sync::Mutex<Option<GpuBuffers>>,
}

/// Persistent GPU work-set, sized for `n` atoms.
struct GpuBuffers {
    n: usize,
    pos_buf: wgpu::Buffer,
    force_buf: wgpu::Buffer,
    aux_buf: wgpu::Buffer,
    staging: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    /// Cell-list (CSR) state: `j_start`/`j_list` buffers + bind group,
    /// allocated on first use and grown as the list grows. `csr_fingerprint`
    /// identifies the list contents the uploaded CSR was built from (rebuild
    /// detection: the CSR is re-uploaded only when the half list changes).
    cell_list: Option<CellListBuffers>,
}

struct CellListBuffers {
    j_start_buf: wgpu::Buffer,
    j_list_buf: wgpu::Buffer,
    j_list_capacity: usize,
    bind_group: wgpu::BindGroup,
    fingerprint: u64,
}

impl GpuPairForceModel {
    /// Acquire a GPU adapter + device and compile the LJ compute shaders.
    ///
    /// Returns `None` (not an error) when no adapter/device is available or a
    /// shader fails to compile — callers should fall back to the CPU model.
    pub fn new(sigma: f64, epsilon: f64, cutoff: f64) -> Option<Self> {
        Self::with_gather(sigma, epsilon, cutoff, GpuGather::AllPairs)
    }

    /// Like [`new`](Self::new), selecting the gather strategy. The GPU-d
    /// ladder (G3) ships [`GpuGather::CellList`] as the performant path and
    /// keeps [`GpuGather::AllPairs`] as the determinism oracle.
    pub fn with_gather(sigma: f64, epsilon: f64, cutoff: f64, gather: GpuGather) -> Option<Self> {
        let sigma = sigma as f32;
        let epsilon = epsilon as f32;
        let cutoff = cutoff as f32;

        // Shifted-force constants at the cutoff (f32), matching the CPU tables.
        let sc = sigma / cutoff;
        let s6c = sc * sc * sc * sc * sc * sc;
        let s12c = s6c * s6c;
        let u_rc = 4.0 * epsilon * (s12c - s6c);
        let f_rc = (24.0 * epsilon / cutoff) * (2.0 * s12c - s6c);

        let params = GpuParams {
            n: 0.0,
            sigma,
            epsilon,
            cutoff,
            cutoff2: cutoff * cutoff,
            u_rc,
            f_rc,
            _pad0: 0.0,
            bx: 0.0,
            by: 0.0,
            bz: 0.0,
            _pad1: 0.0,
        };

        let device_and_queue = pollster::block_on(async {
            let instance = wgpu::Instance::default();
            let adapter = instance
                .request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference: wgpu::PowerPreference::HighPerformance,
                    force_fallback_adapter: false,
                    compatible_surface: None,
                    apply_limit_buckets: false,
                })
                .await
                .ok()?;
            adapter
                .request_device(&wgpu::DeviceDescriptor {
                    label: Some("modena gpu"),
                    required_features: wgpu::Features::empty(),
                    // Native limits (not downlevel): the cell-list kernel binds
                    // 5 storage buffers per stage (positions, forces, aux,
                    // j_start, j_list) vs the downlevel cap of 4. Adapters that
                    // cannot provide this fail the request and the constructor
                    // returns `None` — the same graceful-CPU-fallback path as
                    // a missing adapter.
                    required_limits: wgpu::Limits::default(),
                    experimental_features: Default::default(),
                    memory_hints: Default::default(),
                    trace: wgpu::Trace::default(),
                })
                .await
                .ok()
        })?;
        let (device, queue) = device_and_queue;

        // Compile the shaders, surfacing any WGSL error as `None`.
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("modena lj"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });
        let module_cell_list = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("modena lj csr"),
            source: wgpu::ShaderSource::Wgsl(SHADER_CELL_LIST.into()),
        });

        let uniform_entry = wgpu::BindGroupLayoutEntry {
            binding: 0,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty: wgpu::BufferBindingType::Uniform,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("modena lj bgl"),
            entries: &[
                uniform_entry,
                storage_entry(1, false),
                storage_entry(2, true),
                storage_entry(3, true),
            ],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("modena lj pll"),
            bind_group_layouts: &[Some(&bind_group_layout)],
            immediate_size: 0,
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("modena lj pipeline"),
            layout: Some(&pipeline_layout),
            module: &module,
            entry_point: Some("main"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });

        // Cell-list (CSR) pipeline: same first four bindings plus the two
        // read-only j-range storage buffers.
        let cell_list_bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("modena lj csr bgl"),
            entries: &[
                uniform_entry,
                storage_entry(1, false),
                storage_entry(2, true),
                storage_entry(3, true),
                storage_entry(4, false),
                storage_entry(5, false),
            ],
        });
        let cell_list_pll = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("modena lj csr pll"),
            bind_group_layouts: &[Some(&cell_list_bgl)],
            immediate_size: 0,
        });
        let pipeline_cell_list = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("modena lj csr pipeline"),
            layout: Some(&cell_list_pll),
            module: &module_cell_list,
            entry_point: Some("main"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });

        let uniform = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("modena lj uniform"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        Some(Self {
            device,
            queue,
            pipeline,
            pipeline_cell_list,
            cell_list_bgl,
            uniform,
            params,
            gather,
            buffers: std::sync::Mutex::new(None),
        })
    }

    /// Allocate (or re-allocate, only when `n` changed) the buffers + bind
    /// group. Everything is created once per atom count and reused across
    /// steps — the per-call `create_buffer`×4 + `create_bind_group` churn that
    /// cost ~2 ms/step is gone.
    fn ensure_buffers(&self, n: usize) -> Result<std::sync::MutexGuard<'_, Option<GpuBuffers>>> {
        let mut g = self.buffers.lock().unwrap();
        if g.as_ref().is_none_or(|b| b.n != n) {
            let pos_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("modena gpu pos"),
                size: (3 * n * 4) as wgpu::BufferAddress,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            let force_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("modena gpu forces"),
                size: (3 * n * 4) as wgpu::BufferAddress,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            });
            let aux_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("modena gpu aux"),
                size: (2 * n * 4) as wgpu::BufferAddress,
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                mapped_at_creation: false,
            });
            let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("modena gpu staging"),
                size: (5 * n * 4) as wgpu::BufferAddress,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("modena gpu bg"),
                layout: &self.pipeline.get_bind_group_layout(0),
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: self.uniform.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: pos_buf.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: force_buf.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 3,
                        resource: aux_buf.as_entire_binding(),
                    },
                ],
            });
            *g = Some(GpuBuffers {
                n,
                pos_buf,
                force_buf,
                aux_buf,
                staging,
                bind_group,
                cell_list: None,
            });
        }
        Ok(g)
    }

    /// Convert the half neighbor list into contiguous per-atom j-ranges (CSR):
    /// `j_start[i]..j_start[i+1]` indexes `j_list`, which holds every neighbor
    /// `j` of `i` exactly once (each half-list pair contributes `j` to `i`'s
    /// range and `i` to `j`'s). O(pairs); deterministic order (serial fill).
    fn build_csr(pairs: &[(usize, usize)], n: usize) -> (Vec<u32>, Vec<u32>) {
        assert!(
            2 * pairs.len() < u32::MAX as usize,
            "cell-list j_list exceeds u32 index range"
        );
        let mut j_start = vec![0u32; n + 1];
        for &(i, j) in pairs {
            j_start[i + 1] += 1;
            j_start[j + 1] += 1;
        }
        for i in 0..n {
            j_start[i + 1] += j_start[i];
        }
        let mut j_list = vec![0u32; 2 * pairs.len()];
        let mut cursor = j_start[..n].to_vec();
        for &(i, j) in pairs {
            j_list[cursor[i] as usize] = j as u32;
            cursor[i] += 1;
            j_list[cursor[j] as usize] = i as u32;
            cursor[j] += 1;
        }
        (j_start, j_list)
    }

    /// Cheap identity of the current half list: rebuild count, length, and a
    /// sampled hash of the pair contents. The CSR is re-uploaded only when
    /// this changes (i.e. on neighbor-list rebuilds), not every step.
    fn list_fingerprint(neighbors: &NeighborList, n: usize) -> u64 {
        let pairs = neighbors.pairs();
        let mut h: u64 = 0xcbf2_9ce4_8422_2325; // FNV-1a offset basis
        for x in [
            n,
            pairs.len(),
            neighbors.rebuild_count(),
            neighbors.cutoff().to_bits() as usize,
        ] {
            h ^= x as u64;
            h = h.wrapping_mul(0x1000_0000_01b3);
        }
        // Sample up to 32 pairs spread over the list (first/last included).
        let len = pairs.len();
        let stride = (len / 32).max(1);
        for k in (0..len).step_by(stride) {
            let (i, j) = pairs[k];
            h ^= (i as u64) | ((j as u64) << 32);
            h = h.wrapping_mul(0x1000_0000_01b3);
        }
        h
    }

    /// Make sure the GPU holds a CSR matching the current half list: build +
    /// upload on first use, list change, or capacity growth. Grows the
    /// `j_list` buffer geometrically so steady-state rebuilds only rewrite
    /// the existing allocation.
    fn ensure_cell_list(
        &self,
        buffers: &mut GpuBuffers,
        neighbors: &NeighborList,
        n: usize,
    ) -> Result<()> {
        let fp = Self::list_fingerprint(neighbors, n);
        if let Some(cl) = buffers.cell_list.as_ref() {
            if cl.fingerprint == fp {
                return Ok(());
            }
        }
        let (j_start, j_list) = Self::build_csr(neighbors.pairs(), n);

        let need = j_list.len().max(1); // min 4 bytes: wgpu 0-size storage
        let (j_start_buf, j_list_buf, capacity, bind_group) = match buffers.cell_list.take() {
            Some(cl) if cl.j_list_capacity >= need => {
                self.queue
                    .write_buffer(&cl.j_start_buf, 0, bytemuck::cast_slice(&j_start));
                self.queue
                    .write_buffer(&cl.j_list_buf, 0, bytemuck::cast_slice(&j_list));
                (
                    cl.j_start_buf,
                    cl.j_list_buf,
                    cl.j_list_capacity,
                    cl.bind_group,
                )
            }
            _ => {
                let j_start_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some("modena gpu j_start"),
                    size: (j_start.len() * 4).max(4) as wgpu::BufferAddress,
                    usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                    mapped_at_creation: false,
                });
                let j_list_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
                    label: Some("modena gpu j_list"),
                    size: (need * 4) as wgpu::BufferAddress,
                    usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                    mapped_at_creation: false,
                });
                self.queue
                    .write_buffer(&j_start_buf, 0, bytemuck::cast_slice(&j_start));
                self.queue
                    .write_buffer(&j_list_buf, 0, bytemuck::cast_slice(&j_list));
                let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("modena gpu csr bg"),
                    layout: &self.cell_list_bgl,
                    entries: &[
                        wgpu::BindGroupEntry {
                            binding: 0,
                            resource: self.uniform.as_entire_binding(),
                        },
                        wgpu::BindGroupEntry {
                            binding: 1,
                            resource: buffers.pos_buf.as_entire_binding(),
                        },
                        wgpu::BindGroupEntry {
                            binding: 2,
                            resource: buffers.force_buf.as_entire_binding(),
                        },
                        wgpu::BindGroupEntry {
                            binding: 3,
                            resource: buffers.aux_buf.as_entire_binding(),
                        },
                        wgpu::BindGroupEntry {
                            binding: 4,
                            resource: j_start_buf.as_entire_binding(),
                        },
                        wgpu::BindGroupEntry {
                            binding: 5,
                            resource: j_list_buf.as_entire_binding(),
                        },
                    ],
                });
                (j_start_buf, j_list_buf, need, bind_group)
            }
        };
        buffers.cell_list = Some(CellListBuffers {
            j_start_buf,
            j_list_buf,
            j_list_capacity: capacity,
            bind_group,
            fingerprint: fp,
        });
        Ok(())
    }
}

fn storage_entry(binding: u32, write_only: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage {
                read_only: !write_only,
            },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

impl ForceModel for GpuPairForceModel {
    fn neighbor_cutoff(&self) -> Option<f64> {
        Some(self.params.cutoff as f64)
    }
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut NeighborList,
    ) -> Result<f64> {
        neighbors.validate_pair_scales(system.n())?;
        if neighbors.pair_scales().iter().any(|p| {
            if self.gather == GpuGather::CellList {
                !(p.lj == 1.0 || (p.lj == 0.0 && p.coulomb == 0.0))
            } else {
                p.lj != 1.0
            }
        }) {
            return Err(ModenaError::Config("GPU backend does not support fractional/independently excluded LJ pairs; cell-list gather supports binary 0/0 exclusions".into()));
        }
        if system.box_size.is_triclinic() {
            return Err(ModenaError::Config(
                "GPU pair kernel supports orthorhombic boxes only (triclinic given)".into(),
            ));
        }
        let n = system.n();
        if n == 0 {
            return Ok(0.0);
        }

        // Cell-list gather: the j-ranges come from the caller's Verlet-binned
        // half list (skin-aware rebuild inside `build`), so the CSR upload is
        // amortized like the CPU path. The list must cover the model cutoff.
        if self.gather == GpuGather::CellList {
            neighbors.build(system);
            let list_cutoff = neighbors.cutoff();
            let model_cutoff = self.params.cutoff as f64;
            if list_cutoff + 1e-12 < model_cutoff {
                return Err(ModenaError::Config(format!(
                    "cell-list GPU kernel needs a neighbor list built for at least the model \
                     cutoff {model_cutoff}, got {list_cutoff}"
                )));
            }
        }

        let [bx, by, bz] = system.box_size.size.0;
        let mut params = self.params;
        params.n = n as f32;
        params.bx = bx as f32;
        params.by = by as f32;
        params.bz = bz as f32;
        self.queue
            .write_buffer(&self.uniform, 0, bytemuck::bytes_of(&params));

        // Positions (read-only storage), forces + aux (storage, shader writes).
        let mut pos_f32 = vec![0.0f32; 3 * n];
        for i in 0..n {
            let p = system.pos(i);
            pos_f32[3 * i] = p.0[0] as f32;
            pos_f32[3 * i + 1] = p.0[1] as f32;
            pos_f32[3 * i + 2] = p.0[2] as f32;
        }
        // Persistent buffers + bind group, sized (and only reallocated) for n.
        let mut bg = self.ensure_buffers(n)?;
        let buffers = bg.as_mut().expect("ensure_buffers always leaves Some");
        self.queue
            .write_buffer(&buffers.pos_buf, 0, bytemuck::cast_slice(&pos_f32));

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("modena gpu enc"),
            });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("modena gpu pass"),
                timestamp_writes: None,
            });
            match self.gather {
                GpuGather::AllPairs => {
                    pass.set_pipeline(&self.pipeline);
                    pass.set_bind_group(0, &buffers.bind_group, &[]);
                }
                GpuGather::CellList => {
                    self.ensure_cell_list(buffers, neighbors, n)?;
                    let cl = buffers
                        .cell_list
                        .as_ref()
                        .expect("ensure_cell_list always leaves Some");
                    pass.set_pipeline(&self.pipeline_cell_list);
                    pass.set_bind_group(0, &cl.bind_group, &[]);
                }
            }
            pass.dispatch_workgroups((n as u32).div_ceil(64), 1, 1);
        }
        encoder.copy_buffer_to_buffer(
            &buffers.force_buf,
            0,
            &buffers.staging,
            0,
            (3 * n * 4) as u64,
        );
        encoder.copy_buffer_to_buffer(
            &buffers.aux_buf,
            0,
            &buffers.staging,
            (3 * n * 4) as u64,
            (2 * n * 4) as u64,
        );
        self.queue.submit(std::iter::once(encoder.finish()));

        // Blocking map (mpsc wakes when the copy + map complete).
        let (tx, rx) = std::sync::mpsc::channel();
        buffers
            .staging
            .slice(..)
            .map_async(wgpu::MapMode::Read, move |res| {
                let _ = tx.send(res);
            });
        self.device
            .poll(PollType::wait_indefinitely())
            .map_err(|e| ModenaError::Config(format!("wgpu device poll failed: {e}")))?;
        rx.recv()
            .map_err(|e| ModenaError::Config(format!("wgpu map channel closed: {e}")))?
            .map_err(|e| ModenaError::Config(format!("wgpu buffer map failed: {e}")))?;

        let bytes = buffers.staging.slice(..).get_mapped_range().map_err(|e| {
            // The map succeeded, so release it even when acquiring its view fails.
            buffers.staging.unmap();
            ModenaError::Config(format!("wgpu mapped buffer range failed: {e}"))
        })?;
        let all: &[f32] = bytemuck::cast_slice(&bytes[..]);
        let forces_f32 = &all[..3 * n];
        let aux_f32 = &all[3 * n..5 * n];

        // Accumulate (the trait contract is accumulate, not zero).
        for i in 0..n {
            system.forces[3 * i] += forces_f32[3 * i] as f64;
            system.forces[3 * i + 1] += forces_f32[3 * i + 1] as f64;
            system.forces[3 * i + 2] += forces_f32[3 * i + 2] as f64;
        }
        let mut pe = 0.0f64;
        let mut w = 0.0f64;
        for i in 0..n {
            pe += aux_f32[i] as f64;
            w += aux_f32[n + i] as f64;
        }
        virial.w += w;
        drop(bytes);
        // Required before the next submit reuses the persistent staging
        // buffer (the buffers are now cached across steps).
        buffers.staging.unmap();
        Ok(pe)
    }
}
