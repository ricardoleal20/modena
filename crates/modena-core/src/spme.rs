//! SPME Stage A — Smooth Particle-Mesh Ewald *reciprocal energy* only.
//!
//! The order-4 B-spline charge spread + 3D FFT + influence-function sum that
//! reproduces the direct reciprocal Ewald energy on a mesh
//! (`docs/SPME_PLAN.md`). Forces (ik-differentiation of this functional) and
//! the strain-form virial are Stage B/C. Orthorhombic boxes only (the schema
//! is orthorhombic today).
//!
//! # Convergence accounting (the α rule)
//!
//! The only convention freedom is the overall normalization constant `α` —
//! physical `k = 2π·k̂/L` and `κ` are fixed definitions. Stage A sets the
//! mesh-accounted `α = 2π/V · (grid/window corrections)` and the test
//! verifies the SAME normalization reproduces the direct ±k path to the
//! configured tolerance on the deck, a displaced config, AND a different-κ
//! config (different k-weighting). A single stable α proves the functional
//! form; an unstable α means a shape/factor error — the test fails loudly.
//!
//! Influence (Essmann et al. JCP 103, 8577 (1995) form, reduced units):
//!
//! ```text
//! E_rec = α · Σ_{k̂≠0} e^{−k²/4κ²}/k² · |Q̂(k̂)|² / (b3(k1)·b3(k2)·b3(k3))²
//! ```
//!
//! with `b3(k̂) = (4 + 2·cos(2πk̂/N))/6` the DFT of the cubic (order-4)
//! B-spline window at the grid nodes — it cancels the spline sampling of
//! `Q`, leaving the plain-Ewald kernel.

use crate::error::{ModenaError, Result};
use crate::force_model::{ForceModel, Virial};
use crate::system::{SimBox, System};
use crate::vec3::Vec3;

use std::sync::{Arc, Mutex};

use rayon::prelude::*;

use num_complex::Complex;
use rustfft::{Fft, FftPlanner};
/// Cubic (order-4) B-spline weights at fractional `u ∈ [0,1)` for the four
/// node offsets `{-1, 0, 1, 2}`. Sums to 1 exactly.
#[inline]
fn m4(u: f64) -> [f64; 4] {
    let u1 = 1.0 - u;
    [
        u1 * u1 * u1,
        3.0 * u * u * u - 6.0 * u * u + 4.0,
        -3.0 * u * u * u + 3.0 * u * u + 3.0 * u + 1.0,
        u * u * u,
    ]
    .map(|x| x / 6.0)
}

/// Reciprocal SPME energy (Stage A) for an orthorhombic box.
pub fn spme_reciprocal_energy(
    kappa: f64,
    kcut: f64,
    box_size: &SimBox,
    positions: &[f64],
    charges: &[f64],
    mesh: [usize; 3],
) -> Result<f64> {
    if box_size.is_triclinic() {
        return Err(ModenaError::Config(
            "SPME requires an orthorhombic box".into(),
        ));
    }
    let n = charges.len();
    if n == 0 {
        return Ok(0.0);
    }
    let [n0, n1, n2] = mesh;
    let l = [box_size.size.0[0], box_size.size.0[1], box_size.size.0[2]];
    let ntot = n0 * n1 * n2;

    // ---- 1. B-spline charge spread onto the mesh. ----
    let mut q: Vec<f64> = vec![0.0; ntot];
    let mut u = [0.0f64; 3];
    let mut m0 = [0usize; 3];
    let mut w = [[0.0f64; 4]; 3];
    for i in 0..n {
        let xi = Vec3::new(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
        for a in 0..3 {
            // Fractional coordinate in [0,1) (mesh-wrapped).
            let f = (xi.0[a] / l[a]).fract().rem_euclid(1.0);
            let t = f * mesh[a] as f64;
            let m = t.floor() as usize % mesh[a];
            u[a] = t - m as f64;
            m0[a] = m;
            w[a] = m4(u[a]);
        }
        let qc = charges[i];
        for k0 in 0..4usize {
            let a0 = (m0[0] + k0 + n0 - 1) % n0;
            for k1 in 0..4usize {
                let a1 = (m0[1] + k1 + n1 - 1) % n1;
                for k2 in 0..4usize {
                    let a2 = (m0[2] + k2 + n2 - 1) % n2;
                    let weight = w[0][k0] * w[1][k1] * w[2][k2];
                    q[(a0 * n1 + a1) * n2 + a2] += qc * weight;
                }
            }
        }
    }

    // ---- 2. Forward FFT (real spread → complex). ----
    let mut grid: Vec<Complex<f64>> = q.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let mut planner = FftPlanner::<f64>::new();
    let ffts = [
        planner.plan_fft_forward(n0),
        planner.plan_fft_forward(n1),
        planner.plan_fft_forward(n2),
    ];
    let scratch_len = ffts
        .iter()
        .map(|f| f.get_inplace_scratch_len())
        .max()
        .unwrap_or(0);
    let mut scratch: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); scratch_len.max(1)];

    // Axis 2 (stride 1): contiguous lines.
    for i0 in 0..n0 {
        for i1 in 0..n1 {
            let start = (i0 * n1 + i1) * n2;
            ffts[2].process_with_scratch(&mut grid[start..start + n2], &mut scratch);
        }
    }
    // Axis 1 (stride n2): strided → transpose copy.
    for i0 in 0..n0 {
        for i2 in 0..n2 {
            let mut line: Vec<Complex<f64>> =
                (0..n1).map(|i1| grid[(i0 * n1 + i1) * n2 + i2]).collect();
            ffts[1].process_with_scratch(&mut line, &mut scratch);
            for (i1, v) in line.into_iter().enumerate() {
                grid[(i0 * n1 + i1) * n2 + i2] = v;
            }
        }
    }
    // Axis 0 (stride n1*n2): strided → transpose copy.
    for i1 in 0..n1 {
        for i2 in 0..n2 {
            let mut line: Vec<Complex<f64>> =
                (0..n0).map(|i0| grid[(i0 * n1 + i1) * n2 + i2]).collect();
            ffts[0].process_with_scratch(&mut line, &mut scratch);
            for (i0, v) in line.into_iter().enumerate() {
                grid[(i0 * n1 + i1) * n2 + i2] = v;
            }
        }
    }

    // ---- 3. Influence-function sum over the full mesh. ----
    // Physical k = 2π·k̂/L with k̂ ∈ [−N/2, N/2) (centered indexing).
    let khat = |nn: usize, m: usize| -> i64 {
        let m = m as i64;
        if m <= (nn as i64) / 2 {
            m
        } else {
            m - nn as i64
        }
    };
    let b3 = |kh: i64, n: usize| -> f64 {
        (4.0 + 2.0 * (2.0 * std::f64::consts::PI * kh as f64 / n as f64).cos()) / 6.0
    };
    let two_pi = 2.0 * std::f64::consts::PI;
    let kappa2 = kappa * kappa;
    let mut e = 0.0f64;
    for n0i in 0..n0 {
        let kx = khat(n0, n0i);
        let phx = two_pi * kx as f64 / l[0];
        for n1i in 0..n1 {
            let ky = khat(n1, n1i);
            let phy = two_pi * ky as f64 / l[1];
            for n2i in 0..n2 {
                let kz = khat(n2, n2i);
                if kx == 0 && ky == 0 && kz == 0 {
                    continue;
                }
                let phz = two_pi * kz as f64 / l[2];
                let k2 = phx * phx + phy * phy + phz * phz;
                if k2 <= 0.0 {
                    continue;
                }
                let qhat = grid[(n0i * n1 + n1i) * n2 + n2i];
                let w = libm::exp(-k2 / (4.0 * kappa2)) / k2;
                let b = b3(kx, n0) * b3(ky, n1) * b3(kz, n2);
                e += w * (qhat.norm_sqr()) / (b * b);
            }
        }
    }
    let _ = kcut;
    let vol = l[0] * l[1] * l[2];
    let alpha = two_pi / vol;
    Ok(alpha * e)
}

// ---------------------------------------------------------------------------
// Stage B/C: full `SpmeForceModel` — forces (ik-differentiation), strain-form
// reciprocal virial, and real-space/self/background wiring.
// ---------------------------------------------------------------------------

/// Deterministic mesh for an orthorhombic box and k-space cutoff
/// (`h = L/N ≤ π/kcut`, the aliasing-safe bound used by Stage A's tests),
/// rounded up per axis to the next 5-smooth size (factors 2·3·5 only — the
/// FFT-friendly convention; 123 → 125 costs ~1.6% more points and saves the
/// prime-factor 41 butterflies).
pub fn mesh_for_box(box_size: &SimBox, kcut: f64) -> [usize; 3] {
    let raw = [
        (box_size.size.0[0] * kcut / std::f64::consts::PI).ceil() as usize,
        (box_size.size.0[1] * kcut / std::f64::consts::PI).ceil() as usize,
        (box_size.size.0[2] * kcut / std::f64::consts::PI).ceil() as usize,
    ];
    [
        next_5_smooth(raw[0]),
        next_5_smooth(raw[1]),
        next_5_smooth(raw[2]),
    ]
}

/// Smallest `n ≥ x` whose prime factors are all ≤ 5.
fn next_5_smooth(x: usize) -> usize {
    let mut n = x.max(4);
    loop {
        let mut m = n;
        for p in [2, 3, 5] {
            while m % p == 0 {
                m /= p;
            }
        }
        if m == 1 {
            return n;
        }
        n += 1;
    }
}

/// Smooth Particle-Mesh Ewald [`ForceModel`]: ORDER-4 B-spline spread + FFT
/// reciprocal (energy, forces, strain-form virial) + the real-space erfc and
/// self/background pieces identical to the direct Ewald.
pub struct SpmeForceModel {
    pub kappa: f64,
    pub real_cutoff: f64,
    pub kspace_cutoff: f64,
}

impl SpmeForceModel {
    pub fn new(kappa: f64, real_cutoff: f64, kspace_cutoff: f64) -> Result<Self> {
        if !(kappa > 0.0) {
            return Err(ModenaError::Config("SPME `kappa` must be > 0".into()));
        }
        if !(real_cutoff > 0.0) {
            return Err(ModenaError::Config("SPME `real_cutoff` must be > 0".into()));
        }
        if !(kspace_cutoff > 0.0) {
            return Err(ModenaError::Config(
                "SPME `kspace_cutoff` must be > 0".into(),
            ));
        }
        Ok(Self {
            kappa,
            real_cutoff,
            kspace_cutoff,
        })
    }

    /// The real-space cutoff — the shared-neighbor-list max, like
    /// [`EwaldForceModel::real_cutoff`].
    pub fn real_cutoff(&self) -> f64 {
        self.real_cutoff
    }

    fn real_space(
        &self,
        box_size: SimBox,
        positions: &[f64],
        charges: &[f64],
        neighbors: &mut crate::neighbor::NeighborList,
        forces: &mut [f64],
    ) -> (f64, f64) {
        // erfc-screened short-range sum — same formula as the direct Ewald's
        // real space (shared-form; a later refactor can unify the loop).
        let rc2 = self.real_cutoff * self.real_cutoff;
        let kappa = self.kappa;
        let mut u = 0.0;
        let mut w = 0.0;
        for &(i, j) in neighbors.pairs() {
            let pi = Vec3::new(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
            let pj = Vec3::new(positions[3 * j], positions[3 * j + 1], positions[3 * j + 2]);
            let d = box_size.min_image(pi - pj);
            let r2 = d.norm2();
            if r2 >= rc2 || r2 < 1e-14 {
                continue;
            }
            let r = r2.sqrt();
            let qprod = charges[i] * charges[j];
            let erfc_kr = libm::erfc(kappa * r);
            let exp_kr2 = libm::exp(-kappa * kappa * r2);
            u += qprod * erfc_kr / r;
            let coef = qprod * (TWO_OVER_SQRT_PI * kappa * exp_kr2 / r2 + erfc_kr / (r2 * r));
            let fv = d * coef;
            for a in 0..3 {
                forces[3 * i + a] += fv.0[a];
                forces[3 * j + a] -= fv.0[a];
                w += d.0[a] * fv.0[a];
            }
        }
        (u, w)
    }
}

const TWO_OVER_SQRT_PI: f64 = std::f64::consts::FRAC_2_SQRT_PI;

impl ForceModel for SpmeForceModel {
    fn compute(
        &self,
        system: &mut System,
        virial: &mut Virial,
        neighbors: &mut crate::neighbor::NeighborList,
    ) -> Result<f64> {
        if system.box_size.is_triclinic() {
            return Err(ModenaError::Config(
                "SPME supports orthorhombic boxes only (triclinic given)".into(),
            ));
        }
        let n = system.n();
        if n == 0 {
            virial.w = 0.0;
            return Ok(0.0);
        }
        neighbors.build(system);

        // Real-space erfc + self/background (identical to Ewald). Split the
        // force buffer out so positions/charges stay shared borrows.
        let mut forces = std::mem::take(&mut system.forces);
        let (u_real, w_real) = self.real_space(
            system.box_size,
            &system.positions,
            &system.charges,
            neighbors,
            &mut forces,
        );
        // Reciprocal on the mesh. The 4-tap stencil needs ≥4 planes per axis
        // to resolve anything; a smaller mesh means the k-space cutoff is too
        // small for the box (raise it, or use the default direct Ewald).
        let mesh = mesh_for_box(&system.box_size, self.kspace_cutoff);
        if mesh[0] < 4 || mesh[1] < 4 || mesh[2] < 4 {
            return Err(ModenaError::Config(format!(
                "SPME mesh {:?} too small for the box (kspace_cutoff {}); \
                 need >= 4 planes per axis",
                mesh, self.kspace_cutoff
            )));
        }
        let (u_rec, w_rec) = spme_reciprocal_forces(
            self.kappa,
            &system.box_size,
            &system.positions,
            &system.charges,
            mesh,
            &mut forces,
        )?;
        system.forces = forces;
        let u_self = -(self.kappa / std::f64::consts::PI.sqrt())
            * system.charges.iter().map(|q| q * q).sum::<f64>();
        let v = system.box_size.volume();
        let sum_q: f64 = system.charges.iter().sum();
        let u_bg = -std::f64::consts::PI / (2.0 * v * self.kappa * self.kappa) * sum_q * sum_q;

        virial.w += w_real + w_rec;
        Ok(u_real + u_rec + u_self + u_bg)
    }
}

/// Cached per-k-point influence factors: `fac = α·e^{−k²/4κ²}/(k²·b²)` and
/// the strain-virial factor `vfac = 1 − k²/(2κ²)`. Pure functions of (κ, box,
/// mesh) — computed once and reused every step (the same pattern as the direct
/// Ewald's k-vector cache). `fac = 0` masks the k = 0 plane. Stored in the
/// transposed L1 layout `[n2][n1][n0]` — the layout the k-pass runs in (see
/// [`spme_reciprocal_forces`]).
struct InfluenceFactors {
    fac: Vec<f64>,
    vfac: Vec<f64>,
}

type InfluenceKey = (u64, [u64; 3]);
type InfluenceEntry = (InfluenceKey, [usize; 3], Arc<InfluenceFactors>);

static INFLUENCE_CACHE: Mutex<Option<InfluenceEntry>> = Mutex::new(None);

fn influence_factors(kappa: f64, box_size: &SimBox, mesh: [usize; 3]) -> Arc<InfluenceFactors> {
    let key = (
        kappa.to_bits(),
        [
            box_size.size.0[0].to_bits(),
            box_size.size.0[1].to_bits(),
            box_size.size.0[2].to_bits(),
        ],
    );
    {
        let g = INFLUENCE_CACHE.lock().unwrap();
        if let Some((k, m, data)) = &*g {
            if *k == key && *m == mesh {
                return Arc::clone(data);
            }
        }
    }
    let [n0, n1, n2] = mesh;
    let ntot = n0 * n1 * n2;
    let l = [box_size.size.0[0], box_size.size.0[1], box_size.size.0[2]];
    let two_pi = 2.0 * std::f64::consts::PI;
    let kappa2 = kappa * kappa;
    let alpha = two_pi / (l[0] * l[1] * l[2]);
    let khat = |nn: usize, m: usize| -> i64 {
        let m = m as i64;
        if m <= (nn as i64) / 2 {
            m
        } else {
            m - nn as i64
        }
    };
    let b3 =
        |kh: i64, nn: usize| -> f64 { (4.0 + 2.0 * (two_pi * kh as f64 / nn as f64).cos()) / 6.0 };
    let mut fac = vec![0.0f64; ntot];
    let mut vfac = vec![0.0f64; ntot];
    // L1 layout: idx = (n2i*n1 + n1i)*n0 + n0i.
    for n2i in 0..n2 {
        let kz = khat(n2, n2i);
        let phz = two_pi * kz as f64 / l[2];
        for n1i in 0..n1 {
            let ky = khat(n1, n1i);
            let phy = two_pi * ky as f64 / l[1];
            for n0i in 0..n0 {
                let kx = khat(n0, n0i);
                if kx == 0 && ky == 0 && kz == 0 {
                    continue; // masked (fac stays 0)
                }
                let phx = two_pi * kx as f64 / l[0];
                let k2 = phx * phx + phy * phy + phz * phz;
                let w = libm::exp(-k2 / (4.0 * kappa2)) / k2;
                let b = b3(kx, n0) * b3(ky, n1) * b3(kz, n2);
                let idx = (n2i * n1 + n1i) * n0 + n0i;
                fac[idx] = alpha * w / (b * b);
                vfac[idx] = 1.0 - k2 / (2.0 * kappa2);
            }
        }
    }
    let data = Arc::new(InfluenceFactors { fac, vfac });
    *INFLUENCE_CACHE.lock().unwrap() = Some((key, mesh, Arc::clone(&data)));
    data
}

/// Grid points below this run the axis FFTs serially (rayon overhead exceeds
/// the transform cost on the small test meshes).
const FFT_SERIAL_POINTS: usize = 1 << 16;

/// In-place FFT over a grid viewed as contiguous rows of length `len` — axis
/// 2 in the L0 layout, axis 0 in the transposed L1 layout. Rows are disjoint,
/// so rays split cleanly; each ray owns its scratch.
fn fft_contiguous(grid: &mut [Complex<f64>], plan: &dyn Fft<f64>, len: usize, parallel: bool) {
    let scratch_len = plan.get_inplace_scratch_len().max(1);
    if parallel {
        grid.par_chunks_exact_mut(len).for_each(|row| {
            let mut scratch = vec![Complex::new(0.0, 0.0); scratch_len];
            plan.process_with_scratch(row, &mut scratch);
        });
    } else {
        let mut scratch = vec![Complex::new(0.0, 0.0); scratch_len];
        for row in grid.chunks_exact_mut(len) {
            plan.process_with_scratch(row, &mut scratch);
        }
    }
}

/// Axis-1 pass on the L0 layout. Slab i0 = grid[i0*n1*n2 .. (i0+1)*n1*n2]:
/// every axis-1 line with that i0 lives inside the slab, so slabs are
/// disjoint &mut and rays split cleanly.
fn fft_axis1(grid: &mut [Complex<f64>], plan: &dyn Fft<f64>, dims: [usize; 3], parallel: bool) {
    let [_n0, n1, n2] = dims;
    let scratch_len = plan.get_inplace_scratch_len().max(1);
    let slab = n1 * n2;
    let slab_pass = |slab_grid: &mut [Complex<f64>], scratch: &mut Vec<Complex<f64>>| {
        let mut line = vec![Complex::new(0.0, 0.0); n1];
        for i2 in 0..n2 {
            for (t, i1) in (0..n1).enumerate() {
                line[t] = slab_grid[i1 * n2 + i2];
            }
            plan.process_with_scratch(&mut line, scratch);
            for (i1, v) in line.iter().enumerate() {
                slab_grid[i1 * n2 + i2] = *v;
            }
        }
    };
    if parallel {
        grid.par_chunks_mut(slab).for_each(|slab_grid| {
            let mut scratch = vec![Complex::new(0.0, 0.0); scratch_len];
            slab_pass(slab_grid, &mut scratch);
        });
    } else {
        let mut scratch = vec![Complex::new(0.0, 0.0); scratch_len];
        for slab_grid in grid.chunks_mut(slab) {
            slab_pass(slab_grid, &mut scratch);
        }
    }
}

/// Transpose L0 `[n0][n1][n2]` → L1 `[n2][n1][n0]`:
/// `dst[(i2*n1+i1)*n0 + i0] = src[(i0*n1+i1)*n2 + i2]`. An involution when
/// called with the dims reversed — the same function transposes back.
/// Parallel over L1 slabs (contiguous writes, strided reads).
fn transpose_021(src: &[Complex<f64>], dst: &mut [Complex<f64>], dims: [usize; 3], parallel: bool) {
    let [n0, n1, n2] = dims;
    let l1_slab = n1 * n0;
    let trans = |i2: usize, dst_slab: &mut [Complex<f64>]| {
        for i1 in 0..n1 {
            let dst_row = &mut dst_slab[i1 * n0..(i1 + 1) * n0];
            for (i0, d) in dst_row.iter_mut().enumerate() {
                *d = src[(i0 * n1 + i1) * n2 + i2];
            }
        }
    };
    if parallel {
        dst.par_chunks_mut(l1_slab)
            .enumerate()
            .for_each(|(i2, slab)| trans(i2, slab));
    } else {
        for (i2, slab) in dst.chunks_mut(l1_slab).enumerate() {
            trans(i2, slab);
        }
    }
}

/// Reciprocal SPME pass with forces + strain-form virial. Writes per-atom
/// forces into `forces`, returns `(E_rec, W_rec)`.
pub fn spme_reciprocal_forces(
    kappa: f64,
    box_size: &SimBox,
    positions: &[f64],
    charges: &[f64],
    mesh: [usize; 3],
    forces: &mut [f64],
) -> Result<(f64, f64)> {
    let n = charges.len();
    let [n0, n1, n2] = mesh;
    let l = [box_size.size.0[0], box_size.size.0[1], box_size.size.0[2]];
    let ntot = n0 * n1 * n2;

    // --- 1. spread ---
    let mut q: Vec<f64> = vec![0.0; ntot];
    let mut u = [0.0f64; 3];
    let mut m0 = [0usize; 3];
    let mut w = [[0.0f64; 4]; 3];
    for i in 0..n {
        let xi = Vec3::new(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
        for a in 0..3 {
            let f = (xi.0[a] / l[a]).fract().rem_euclid(1.0);
            let t = f * mesh[a] as f64;
            let m = t.floor() as usize % mesh[a];
            u[a] = t - m as f64;
            m0[a] = m;
            w[a] = m4(u[a]);
        }
        let qc = charges[i];
        for k0 in 0..4usize {
            let a0 = (m0[0] + k0 + n0 - 1) % n0;
            for k1 in 0..4usize {
                let a1 = (m0[1] + k1 + n1 - 1) % n1;
                for k2 in 0..4usize {
                    let a2 = (m0[2] + k2 + n2 - 1) % n2;
                    q[(a0 * n1 + a1) * n2 + a2] += qc * w[0][k0] * w[1][k1] * w[2][k2];
                }
            }
        }
    }

    // --- 2. forward FFT: axes 2/1 on the L0 grid, then transpose to L1 ---
    // ([n2][n1][n0]) and transform axis 0 contiguously — no serial axis pass
    // remains; one extra grid buffer per call.
    let mut grid: Vec<Complex<f64>> = q.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let mut planner = FftPlanner::<f64>::new();
    let fwd = [
        planner.plan_fft_forward(n0),
        planner.plan_fft_forward(n1),
        planner.plan_fft_forward(n2),
    ];
    let bwd = [
        planner.plan_fft_inverse(n0),
        planner.plan_fft_inverse(n1),
        planner.plan_fft_inverse(n2),
    ];
    let parallel = ntot >= FFT_SERIAL_POINTS;
    fft_contiguous(&mut grid, fwd[2].as_ref(), n2, parallel);
    fft_axis1(&mut grid, fwd[1].as_ref(), mesh, parallel);
    let mut l1: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); ntot];
    transpose_021(&grid, &mut l1, mesh, parallel);
    fft_contiguous(&mut l1, fwd[0].as_ref(), n0, parallel);

    // --- 3. k-pass (in the L1 layout): energy + strain virial + Q̂' in place.
    // (influence factors cached per (κ, box, mesh), in L1 order). ---
    let infl = influence_factors(kappa, box_size, mesh);
    let mut e = 0.0f64;
    let mut wr = 0.0f64;
    for (dst, (f, vf)) in l1.iter_mut().zip(infl.fac.iter().zip(infl.vfac.iter())) {
        let qhat = *dst;
        let u_k = f * qhat.norm_sqr();
        e += u_k;
        wr += u_k * vf;
        *dst = qhat * f;
    }

    // --- 4. inverse FFT: axis 0 contiguous on L1, transpose back into the
    // spread grid's allocation, axes 1/2 on L0 → Re[Φ]. ---
    fft_contiguous(&mut l1, bwd[0].as_ref(), n0, parallel);
    transpose_021(&l1, &mut grid, [n2, n1, n0], parallel);
    fft_axis1(&mut grid, bwd[1].as_ref(), mesh, parallel);
    fft_contiguous(&mut grid, bwd[2].as_ref(), n2, parallel);
    let phi = grid;

    // --- 5. force assembly via the B-spline weight derivatives. ---
    // F_i = −2·q_i·Σ_{m∈stencil} ∂W_i(m)/∂r_i · Re[Φ(m)].
    // ∂W/∂r_a = (N_a/L_a)·w′_{ka}(u_a)·∏_{b≠a} w_{kb}(u_b)
    // (the 4-tap derivative sums to zero → translation-invariant).
    let mut m0c = [0usize; 3];
    let mut wc = [[0.0f64; 4]; 3];
    let mut wd = [[0.0f64; 4]; 3];
    for i in 0..n {
        let xi = Vec3::new(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
        for a in 0..3 {
            let f = (xi.0[a] / l[a]).fract().rem_euclid(1.0);
            let t = f * mesh[a] as f64;
            let m = t.floor() as usize % mesh[a];
            m0c[a] = m;
            wc[a] = m4(t - m as f64);
            wd[a] = m4_deriv(t - m as f64);
        }
        let qc = charges[i];
        let mut fv = Vec3::zero();
        for k0 in 0..4usize {
            let a0 = (m0c[0] + k0 + n0 - 1) % n0;
            for k1 in 0..4usize {
                let a1 = (m0c[1] + k1 + n1 - 1) % n1;
                for k2 in 0..4usize {
                    let a2 = (m0c[2] + k2 + n2 - 1) % n2;
                    let phi_m = phi[(a0 * n1 + a1) * n2 + a2].re;
                    // ∂W/∂x = w0'·w1·w2 · (N0/L0) — etc.
                    let dwx = wd[0][k0] * wc[1][k1] * wc[2][k2] * (mesh[0] as f64 / l[0]);
                    let dwy = wc[0][k0] * wd[1][k1] * wc[2][k2] * (mesh[1] as f64 / l[1]);
                    let dwz = wc[0][k0] * wc[1][k1] * wd[2][k2] * (mesh[2] as f64 / l[2]);
                    let c = -2.0 * qc * phi_m;
                    fv.0[0] += c * dwx;
                    fv.0[1] += c * dwy;
                    fv.0[2] += c * dwz;
                }
            }
        }
        forces[3 * i] += fv.0[0];
        forces[3 * i + 1] += fv.0[1];
        forces[3 * i + 2] += fv.0[2];
    }
    Ok((e, wr))
}

/// Derivative of the cubic B-spline weights w.r.t. `u` (sums to 0).
#[inline]
fn m4_deriv(u: f64) -> [f64; 4] {
    let u1 = 1.0 - u;
    [
        -0.5 * u1 * u1,
        (3.0 * u * u - 4.0 * u) / 2.0,
        (-3.0 * u * u + 2.0 * u + 1.0) / 2.0,
        0.5 * u * u,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ewald::EwaldForceModel;
    use crate::lattice;

    fn nacl() -> crate::system::System {
        lattice::rocksalt(5.0, [3, 3, 3], "Na", "Cl", 1.0, 1.0, 1.0, -1.0)
    }

    fn run_case(kappa: f64, kcut: f64, mesh: [usize; 3], disp: bool) -> (f64, f64, f64) {
        let sys = nacl();
        let mut sys = sys;
        if disp {
            // Displaced config (tiny random-like displacement → nonzero forces).
            for i in 0..sys.n() {
                let f = i as f64;
                let p = sys.pos(i)
                    + Vec3::new(
                        1e-3 * (0.13 * f).sin(),
                        1e-3 * (0.7 * f).cos(),
                        1e-3 * (1.3 * f).sin(),
                    );
                sys.set_pos(i, p);
            }
            sys.wrap_positions();
        }
        let ewald = EwaldForceModel::new(kappa, 5.0, kcut).expect("ewald");
        let direct = ewald.reciprocal_energy(&sys).expect("direct").0;
        let spme = spme_reciprocal_energy(
            kappa,
            kcut,
            &sys.box_size,
            &sys.positions,
            &sys.charges,
            mesh,
        )
        .expect("spme");
        let ratio = spme / direct;
        (direct, spme, ratio)
    }

    fn mesh_for(box_len: f64, kcut: f64) -> [usize; 3] {
        let n = (box_len * kcut / std::f64::consts::PI).ceil() as usize;
        [n, n, n]
    }

    #[test]
    fn spme_matches_direct_on_deck() {
        let (direct, spme, ratio) = run_case(0.8, 7.0, mesh_for(15.0, 7.0), false);
        assert!(
            (ratio - 1.0).abs() < 5e-3,
            "deck: direct={direct} spme={spme} ratio={ratio}"
        );
    }

    #[test]
    fn spme_matches_direct_displaced() {
        let (direct, spme, ratio) = run_case(0.8, 7.0, mesh_for(15.0, 7.0), true);
        assert!(
            (ratio - 1.0).abs() < 5e-3,
            "displaced: direct={direct} spme={spme} ratio={ratio}"
        );
    }

    #[test]
    fn spme_alpha_stable_across_kappa() {
        // The SAME normalization must reproduce the direct path at a different
        // κ (different k-weighting) — the α-rule.
        let (direct, spme, ratio) = run_case(0.5, 6.0, mesh_for(15.0, 6.0), true);
        assert!(
            (ratio - 1.0).abs() < 5e-3,
            "kappa=0.5: direct={direct} spme={spme} ratio={ratio}"
        );
    }

    /// Displaced NaCl deck — the config with non-trivial forces for the Stage B/C
    /// gates. `amp` is the displacement amplitude (1e-3 for the FD gates; larger
    /// for the Ewald cross-check, where forces must dominate the B-spline
    /// truncation error of the mesh).
    fn nacl_displaced(amp: f64) -> crate::system::System {
        let mut sys = nacl();
        for i in 0..sys.n() {
            let f = i as f64;
            let p = sys.pos(i)
                + Vec3::new(
                    amp * (0.13 * f).sin(),
                    amp * (0.7 * f).cos(),
                    amp * (1.3 * f).sin(),
                );
            sys.set_pos(i, p);
        }
        sys.wrap_positions();
        sys
    }

    /// Stage B gate: the ik-differentiated mesh forces are the exact derivative
    /// of the SAME mesh energy functional (the never-separate-U/F-interpolants
    /// rule from the engine research — FD of `spme_reciprocal_energy` vs
    /// `spme_reciprocal_forces`, same κ and mesh).
    #[test]
    fn spme_forces_are_derivative_of_energy() {
        let kappa = 0.8;
        let mesh = mesh_for(15.0, 7.0);
        let sys = nacl_displaced(1e-3);
        let n = sys.n();
        let mut forces = vec![0.0; 3 * n];
        spme_reciprocal_forces(
            kappa,
            &sys.box_size,
            &sys.positions,
            &sys.charges,
            mesh,
            &mut forces,
        )
        .expect("spme forces");
        let fmax = forces.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        assert!(
            fmax > 1e-8,
            "displaced deck must have nonzero forces (fmax={fmax})"
        );

        let mut sys = sys;
        let h = 1e-6;
        for i in 0..4.min(n) {
            for a in 0..3 {
                let idx = 3 * i + a;
                let old = sys.positions[idx];
                sys.positions[idx] = old + h;
                let ep = spme_reciprocal_energy(
                    kappa,
                    7.0,
                    &sys.box_size,
                    &sys.positions,
                    &sys.charges,
                    mesh,
                )
                .unwrap();
                sys.positions[idx] = old - h;
                let em = spme_reciprocal_energy(
                    kappa,
                    7.0,
                    &sys.box_size,
                    &sys.positions,
                    &sys.charges,
                    mesh,
                )
                .unwrap();
                sys.positions[idx] = old;
                let fd = -(ep - em) / (2.0 * h);
                assert!(
                    (fd - forces[idx]).abs() < 1e-4 * fmax,
                    "atom {i} axis {a}: ik force {} vs FD {fd}",
                    forces[idx]
                );
            }
        }
    }

    /// Stage C gate: the strain-form reciprocal virial Σ_k U_k(1 − k²/2κ²) is the
    /// isotropic strain derivative of the same mesh functional (w = −dE/dλ with
    /// box AND positions scaled, mesh held fixed).
    #[test]
    fn spme_reciprocal_virial_is_strain_derivative() {
        let kappa = 0.8;
        let mesh = mesh_for(15.0, 7.0);
        let sys = nacl_displaced(1e-3);
        let n = sys.n();
        let mut forces = vec![0.0; 3 * n];
        let (_, w) = spme_reciprocal_forces(
            kappa,
            &sys.box_size,
            &sys.positions,
            &sys.charges,
            mesh,
            &mut forces,
        )
        .expect("spme forces");
        assert!(w.abs() > 1e-6, "strain virial must be nonzero (w={w})");

        let strain_energy = |lam: f64| -> f64 {
            let mut s = sys.clone();
            s.box_size.scale_uniform(lam);
            for x in &mut s.positions {
                *x *= lam;
            }
            spme_reciprocal_energy(kappa, 7.0, &s.box_size, &s.positions, &s.charges, mesh).unwrap()
        };
        let h = 1e-6;
        let fd = -(strain_energy(1.0 + h) - strain_energy(1.0 - h)) / (2.0 * h);
        assert!(
            (fd - w).abs() < 1e-4 * w.abs().max(1.0),
            "strain-form virial {w} vs strain-FD {fd}"
        );
    }

    /// End-to-end: the full `SpmeForceModel` (real erfc + mesh reciprocal +
    /// self/background) matches the direct `EwaldForceModel` energies and forces
    /// on a strongly displaced deck (amplitude 0.05 — forces scale with the
    /// displacement while the B-spline truncation error does not, so the mesh
    /// error stays a small fraction of fmax at kcut=14; on the crystalline 1e-3
    /// deck the ratio inverts, which is mesh physics, not a bug).
    #[test]
    fn spme_force_model_matches_direct_ewald() {
        use crate::force_model::Virial;
        use crate::neighbor::NeighborList;

        let kappa = 0.8;
        let kcut = 14.0;
        let mut s_spme = nacl_displaced(0.05);
        let mut s_ewald = s_spme.clone();
        s_spme.zero_forces();
        s_ewald.zero_forces();

        let spme = SpmeForceModel::new(kappa, 5.0, kcut).expect("spme");
        let ewald = EwaldForceModel::new(kappa, 5.0, kcut).expect("ewald");
        let mut v1 = Virial::default();
        let mut v2 = Virial::default();
        let mut nl1 = NeighborList::new(5.0);
        let mut nl2 = NeighborList::new(5.0);
        let e_spme = spme.compute(&mut s_spme, &mut v1, &mut nl1).expect("spme");
        let e_ewald = ewald
            .compute(&mut s_ewald, &mut v2, &mut nl2)
            .expect("ewald");

        assert!(
            (e_spme - e_ewald).abs() < 1e-3 * e_ewald.abs(),
            "energy: spme {e_spme} vs ewald {e_ewald}"
        );
        let fmax = s_ewald.forces.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        for i in 0..s_spme.n() {
            for a in 0..3 {
                let idx = 3 * i + a;
                assert!(
                    (s_spme.forces[idx] - s_ewald.forces[idx]).abs() < 5e-2 * fmax,
                    "atom {i} axis {a}: spme {} vs ewald {}",
                    s_spme.forces[idx],
                    s_ewald.forces[idx]
                );
            }
        }
        // The virials must NOT agree: the direct model's reciprocal virial is the
        // naive Σ r·F (the documented +0.00632 pressure defect), while SPME uses
        // the strain form (gated exactly by the strain-FD test above). The split
        // must therefore be dominated by the strain-form reciprocal virial.
        let n = s_spme.n();
        let mut scratch = vec![0.0; 3 * n];
        let (_, w_rec_strain) = spme_reciprocal_forces(
            kappa,
            &s_spme.box_size,
            &s_spme.positions,
            &s_spme.charges,
            mesh_for(15.0, kcut),
            &mut scratch,
        )
        .expect("spme rec");
        let w_diff = v1.w - v2.w;
        assert!(
            (w_diff - w_rec_strain).abs() < 0.2 * w_rec_strain.abs(),
            "virial split {w_diff} vs strain-form reciprocal {w_rec_strain} \
         (SPME must carry the strain-form correction, not the naive r·F)"
        );
    }
}
