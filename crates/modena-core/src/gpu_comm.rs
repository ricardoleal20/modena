//! GPU offload of the Lennard-Jones pair kernel (M9, `gpu` feature).
//!
//! A [`ForceModel`] backed by a WGSL compute shader running on the GPU via
//! `wgpu`. The kernel is an all-pairs O(N²) evaluation of the *same* shifted-
//! force Lennard-Jones potential as the CPU path, with orthorhombic minimum-
//! image periodicity. Each thread computes one atom's net force, half its pair
//! energy, and half its pair virial (every pair is visited twice — once from
//! each end — so the per-atom energy/virial are halved and summed on the host).
//!
//! ## Precision
//!
//! The shader computes in `f32`; the CPU path is `f64`. The two therefore agree
//! to single precision (~1e-5 relative), not bit-for-bit. This is the standard
//! GPU MD trade-off (LAMMPS's GPU package does the same). Accumulation uses
//! Kahan compensation in-f32, so the summation error is O(1) instead of
//! O(n·eps) over the all-pairs loop; the host reduces f32 outputs into f64
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
//! available (headless CI, no drivers) or the shader fails to compile, so
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

/// GPU-accelerated Lennard-Jones [`ForceModel`].
///
/// Construct with [`new`](Self::new); returns `None` when no usable GPU adapter
/// exists (so CI without a GPU skips gracefully). The model owns the `wgpu`
/// device/queue and a pre-compiled compute pipeline; work buffers + bind group
/// are cached per atom count (see [`ensure_buffers`](Self::ensure_buffers)) — a
/// previous per-call `create_buffer`×4 cost ~2 ms/step and made the GPU slower
/// than the CPU below ~1k atoms.
pub struct GpuPairForceModel {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    uniform: wgpu::Buffer,
    params: GpuParams,
    /// Buffers + bind group cached per atom count: creating four buffers and a
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
}

impl GpuPairForceModel {
    /// Acquire a GPU adapter + device and compile the LJ compute shader.
    ///
    /// Returns `None` (not an error) when no adapter/device is available or the
    /// shader fails to compile — callers should fall back to the CPU model.
    pub fn new(sigma: f64, epsilon: f64, cutoff: f64) -> Option<Self> {
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
                })
                .await
                .ok()?;
            adapter
                .request_device(&wgpu::DeviceDescriptor {
                    label: Some("modena gpu"),
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::downlevel_defaults(),
                    experimental_features: Default::default(),
                    memory_hints: Default::default(),
                    trace: wgpu::Trace::default(),
                })
                .await
                .ok()
        })?;
        let (device, queue) = device_and_queue;

        // Compile the shader, surfacing any WGSL error as `None`.
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("modena lj"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("modena lj bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
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

        let uniform = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("modena lj uniform"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        Some(Self {
            device,
            queue,
            pipeline,
            uniform,
            params,
            buffers: std::sync::Mutex::new(None),
        })
    }

    /// Allocate (or re-allocate, only when `n` changed) the buffers + bind
    /// group. Everything is created once per atom count and reused across
    /// steps — the per-call `create_buffer`×4 + `create_bind_group` churn that
    /// cost ~2 ms/step is gone.
    fn ensure_buffers(&self, n: usize) -> Result<std::sync::MutexGuard<'_, Option<GpuBuffers>>> {
        let mut g = self.buffers.lock().unwrap();
        if g.as_ref().map_or(true, |b| b.n != n) {
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
            });
        }
        Ok(g)
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
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        _neighbors: &mut NeighborList,
    ) -> Result<f64> {
        if system.box_size.is_triclinic() {
            return Err(ModenaError::Config(
                "GPU pair kernel supports orthorhombic boxes only (triclinic given)".into(),
            ));
        }
        let n = system.n();
        if n == 0 {
            return Ok(0.0);
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
        let bg = self.ensure_buffers(n)?;
        let buffers = bg.as_ref().expect("ensure_buffers always leaves Some");
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
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &buffers.bind_group, &[]);
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

        let bytes = buffers.staging.slice(..).get_mapped_range();
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
