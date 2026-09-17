//! Full-precision measurement harness for the bench/validation decks.
//!
//! Two modes:
//!
//! `ewald-decompose` — reproduces the Ewald validation system (rocksalt
//! a=5.0, 3×3×3, charges ±1, Buckingham A=1000/B=1/rho=0.3 cut 5.0, Ewald
//! κ=0.5/rc=5.0/kc=5.0) and decomposes the potential energy at the exact
//! lattice vs a minimally displaced configuration, comparing the engine's
//! Ewald evaluation against an independent brute-force O(N²) real-space +
//! direct reciprocal-space reference implemented here.
//!
//! `shake-measure` — reproduces the SHAKE validation system (12-atom chain,
//! LJ σ=1/ε=0.1/cut=2.5, 11 rigid bonds of length 1.0) by mirroring the
//! runner's exact construction and stepping loop, measuring per-step
//! bond-length residuals |r−d| and RATTLE velocity residuals |r·v| at full
//! f64 precision over the whole run.

use std::path::Path;

use modena_core::{
    rocksalt, Bond, BondedForceModel, BondedTopology, Buckingham, BuckinghamOverride,
    BuckinghamTypeParams, CompositeForceModel, Constraint, Constraints, EwaldForceModel,
    ForceModel, LennardJones, NeighborList, PairForceModel, SimBox, System, Vec3, VelocityVerlet,
    Virial, Integrator,
};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};

/// Exact replica of `runner::init_velocities` (private there; duplicated here
/// so the harness reproduces the runner's run bit-for-bit).
fn init_velocities(system: &mut System, temperature: f64, seed: u64) {
    let mut rng = StdRng::seed_from_u64(seed);
    let n = system.n();
    let normal = StandardNormal;
    let mut p_com = Vec3::zero();
    let mut total_m = 0.0;
    for i in 0..n {
        let m = system.masses[i];
        let s = (temperature / m).sqrt();
        let gx: f64 = normal.sample(&mut rng);
        let gy: f64 = normal.sample(&mut rng);
        let gz: f64 = normal.sample(&mut rng);
        let v = Vec3::new(s * gx, s * gy, s * gz);
        system.set_vel(i, v);
        p_com += v * m;
        total_m += m;
    }
    let v_com = p_com * (1.0 / total_m);
    for i in 0..n {
        let v = system.vel(i) - v_com;
        system.set_vel(i, v);
    }
    let ke = system.kinetic_energy();
    let dof = system.dof() as f64;
    if ke > 0.0 && dof > 0.0 {
        let current_t = 2.0 * ke / dof;
        if current_t > 0.0 {
            let lambda = (temperature / current_t).sqrt();
            for i in 0..n {
                let v = system.vel(i) * lambda;
                system.set_vel(i, v);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Ewald decomposition
// ---------------------------------------------------------------------------

/// Brute-force real-space Ewald: Σ_{i<j, r<rc} q_i q_j erfc(κr)/r over ALL
/// pairs (no cell list), minimum image. Independent of the engine's neighbor
/// list. Returns (energy, virial Σ r·F).
fn brute_real(system: &System, kappa: f64, rc: f64) -> (f64, f64) {
    let n = system.n();
    let rc2 = rc * rc;
    let mut u = 0.0;
    let mut w = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = system.box_size.min_image(system.pos(i) - system.pos(j));
            let r2 = d.norm2();
            if r2 >= rc2 || r2 < 1e-14 {
                continue;
            }
            let r = r2.sqrt();
            let qprod = system.charge(i) * system.charge(j);
            let e = erfc(kappa * r);
            u += qprod * e / r;
            // F_i = d · qprod [ (2κ/√π) e^{−κ²r²}/r² + erfc(κr)/r³ ]  (engine formula)
            let coef = qprod
                * (2.0 * kappa / std::f64::consts::PI.sqrt() * (-kappa * kappa * r2).exp()
                    / r2
                    + e / (r2 * r));
            w += d.dot(d * coef); // r · F
        }
    }
    (u, w)
}

/// erfc via the standardAbramowitz–Stegun-free route: use the incomplete
/// gamma identity erfc(x) = Γ(1/2, x²)/√π with a Continued-fraction-free
/// series would be slow; instead shell out to libm through a tiny FFI-free
/// numeric route: erfc(x) = exp(-x²)·erfcx-free... Simplest robust approach:
/// use the engine's own libm through modena? Not exported. Use the classic
/// Numerical-Recipes-style rational approximation, accurate to ~1e-12? That
/// is not enough for a 1e-10-level comparison — but our target discrepancy is
/// 2.6e-2, so 1e-12 accuracy is ample. Actually simpler: erfc(x) for x ≥ 2.5
/// (our regime: κr = 0.5·[0,5] → x ∈ [0, 2.5]) — use a high-accuracy series.
///
/// We use the confluent-Hypergeometric-free double-precision erfc:
/// for x < 2.0 the Taylor series of erf, otherwise the asymptotic continued
/// fraction. To keep this short and robust we instead integrate the Gaussian
/// with Gauss–Legendre? Overkill.
///
/// Pragmatic choice: rational approximation (Abramowitz & Stegun 7.1.26 has
/// 1.5e-7 accuracy — too coarse). We need better: use the series
/// erf(x) = 2/√π Σ (-1)^n x^(2n+1)/(n!(2n+1)) for x ≤ 2.5 (converges, terms
/// to n≈30, accuracy ~1e-16 in double for these x), and erfc = 1 − erf.
fn erf_series(x: f64) -> f64 {
    // erf(x) = 2/sqrt(pi) * sum_{n>=0} (-1)^n x^{2n+1} / (n! (2n+1))
    let xx = x * x;
    let mut term = x;
    let mut sum = x;
    for n in 1..80 {
        term *= -xx / n as f64;
        let add = term / (2.0 * n as f64 + 1.0);
        sum += add;
        if add.abs() < 1e-18 * sum.abs() {
            break;
        }
    }
    2.0 / std::f64::consts::PI.sqrt() * sum
}

fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        2.0 - erfc(-x)
    } else {
        1.0 - erf_series(x)
    }
}

/// Brute-force reciprocal-space Ewald:
/// U_rec = (2π/V) Σ_{k≠0, |k|<kc} e^{−k²/(4κ²)}/k² · |Σ_j q_j e^{ik·r_j}|².
/// Also returns the correct strain virial W_rec = Σ_k U_k·(1 − k²/(2κ²))
/// (κ fixed under strain — the LAMMPS convention with user-set gewald).
fn brute_reciprocal(system: &System, kappa: f64, kc: f64) -> (f64, f64, usize) {
    let v = system.box_size.volume();
    let l = [
        system.box_size.size.0[0],
        system.box_size.size.0[1],
        system.box_size.size.0[2],
    ];
    let two_pi = 2.0 * std::f64::consts::PI;
    let n_max = (l[0].max(l[1]).max(l[2]) * kc / two_pi).ceil() as i32;
    let kc2 = kc * kc;
    let kappa2 = kappa * kappa;
    let n = system.n();
    let mut u = 0.0;
    let mut w_strain = 0.0;
    let mut count = 0usize;
    for n0 in -n_max..=n_max {
        for n1 in -n_max..=n_max {
            for n2 in -n_max..=n_max {
                if n0 == 0 && n1 == 0 && n2 == 0 {
                    continue;
                }
                let k = [
                    two_pi * n0 as f64 / l[0],
                    two_pi * n1 as f64 / l[1],
                    two_pi * n2 as f64 / l[2],
                ];
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                if k2 >= kc2 || k2 <= 0.0 {
                    continue;
                }
                count += 1;
                let w = (-k2 / (4.0 * kappa2)).exp() / k2;
                let mut re = 0.0;
                let mut im = 0.0;
                for j in 0..n {
                    let rj = system.pos(j);
                    let phase = k[0] * rj.0[0] + k[1] * rj.0[1] + k[2] * rj.0[2];
                    re += system.charge(j) * phase.cos();
                    im += system.charge(j) * phase.sin();
                }
                let uk = two_pi / v * w * (re * re + im * im);
                u += uk;
                w_strain += uk * (1.0 - k2 / (2.0 * kappa2));
            }
        }
    }
    (u, w_strain, count)
}

fn brute_self(system: &System, kappa: f64) -> f64 {
    let sum_q2: f64 = system.charges.iter().map(|q| q * q).sum();
    -(kappa / std::f64::consts::PI.sqrt()) * sum_q2
}

fn brute_bg(system: &System, kappa: f64) -> f64 {
    let sum_q: f64 = system.charges.iter().sum();
    -std::f64::consts::PI / (2.0 * system.box_size.volume() * kappa * kappa) * sum_q * sum_q
}

/// One line of the Ewald decomposition table.
struct Decomp {
    label: &'static str,
    engine_ewald: f64,
    engine_w_virial: f64,
    energy_only: f64,
    bf_real: f64,
    bf_real_w: f64,
    bf_rec: f64,
    bf_rec_wstrain: f64,
    bf_self: f64,
    bf_bg: f64,
    bf_total: f64,
    n_k: usize,
}

fn decompose(system: &System, ewald: &EwaldForceModel, label: &'static str) -> Decomp {
    // Engine Ewald via the ForceModel compute path (forces + energy + virial).
    let mut sys = system.clone();
    sys.zero_forces();
    let mut virial = Virial::default();
    let mut nl = NeighborList::new(ewald.real_cutoff);
    let engine_ewald = ewald.compute(&mut sys, &mut virial, &mut nl).unwrap();
    let engine_w_virial = virial.w;

    // Engine Ewald via the energy-only path.
    let mut nl2 = NeighborList::new(ewald.real_cutoff);
    let energy_only = ewald.energy(system, &mut nl2).unwrap();

    // Brute-force reference.
    let (bf_real, bf_real_w) = brute_real(system, ewald.kappa, ewald.real_cutoff);
    let (bf_rec, bf_rec_wstrain, n_k) =
        brute_reciprocal(system, ewald.kappa, ewald.kspace_cutoff);
    let bf_self = brute_self(system, ewald.kappa);
    let bf_bg = brute_bg(system, ewald.kappa);

    Decomp {
        label,
        engine_ewald,
        engine_w_virial,
        energy_only,
        bf_real,
        bf_real_w,
        bf_rec,
        bf_rec_wstrain,
        bf_self,
        bf_bg,
        bf_total: bf_real + bf_rec + bf_self + bf_bg,
        n_k,
    }
}

fn ewald_decompose() {
    // Exactly the bench/validation/ewald/modena.toml system.
    let mut system = rocksalt(5.0, [3, 3, 3], "Na", "Cl", 1.0, 1.0, 1.0, -1.0);
    init_velocities(&mut system, 0.1, 87287);

    let ewald = EwaldForceModel::new(0.8, 5.0, 7.0).unwrap();

    // Buckingham pair model (same parameters as the TOML).
    let buck = Buckingham::from_types(
        &[
            BuckinghamTypeParams {
                a: 1000.0,
                b: 1.0,
                rho: 0.3,
            },
            BuckinghamTypeParams {
                a: 1000.0,
                b: 1.0,
                rho: 0.3,
            },
        ],
        modena_core::MixingRule::LorentzBerthelot,
        &[BuckinghamOverride {
            types: [0, 1],
            a: 1000.0,
            b: 1.0,
            rho: 0.3,
        }],
        5.0,
    )
    .unwrap();
    let pair = PairForceModel::new(buck);

    let n = system.n();
    println!("atoms: {n}");

    // --- Configuration A: the exact lattice (step-0 state). ---
    let mut sys0 = system.clone();
    let mut nl = NeighborList::new(5.0);
    let mut virial = Virial::default();
    sys0.zero_forces();
    let pe_pair0 = pair.compute(&mut sys0, &mut virial, &mut nl).unwrap();
    let d0 = decompose(&system, &ewald, "exact lattice (x0)");

    // --- Configuration B: displaced by v·dt with dt = 1e-6 (the kick does
    // not change v at the exact lattice because F = 0 there, so this is
    // exactly the runner's step-1 configuration). ---
    let dt = 1e-6_f64;
    let mut sys1 = system.clone();
    for i in 0..n {
        let x = sys1.pos(i) + sys1.vel(i) * dt;
        sys1.set_pos(i, x);
    }
    sys1.wrap_positions();
    let mut nl1 = NeighborList::new(5.0);
    let mut virial1 = Virial::default();
    sys1.zero_forces();
    let pe_pair1 = pair.compute(&mut sys1, &mut virial1, &mut nl1).unwrap();
    let d1 = decompose(&sys1, &ewald, "displaced x0+v*1e-6 (x1)");

    // --- Report. ---
    for d in [&d0, &d1] {
        println!();
        println!("== {} ==", d.label);
        println!("  engine compute() Ewald PE : {:.12}", d.engine_ewald);
        println!("  engine energy()   Ewald PE : {:.12}", d.energy_only);
        println!("  brute real        : {:.12}", d.bf_real);
        println!("  brute reciprocal  : {:.12}  ({} k-vectors)", d.bf_rec, d.n_k);
        println!("  brute self        : {:.12}", d.bf_self);
        println!("  brute background  : {:.12}", d.bf_bg);
        println!("  brute total       : {:.12}", d.bf_total);
        println!(
            "  engine - brute    : {:.12}",
            d.engine_ewald - d.bf_total
        );
        println!("  engine Ewald virial W (Σ r·F, incl. reciprocal): {:.12}", d.engine_w_virial);
        println!("  brute real-space virial W_real (Σ r·F): {:.12}", d.bf_real_w);
        println!("  engine reciprocal virial (Σ r·F form)  : {:.12}", d.engine_w_virial - d.bf_real_w);
        println!(
            "  reciprocal strain virial W_rec = Σ U_k(1-k²/2κ²): {:.12}",
            d.bf_rec_wstrain
        );
    }

    println!();
    println!("== Deltas (x1 - x0) ==");
    println!("  pair PE       : {:.12}", pe_pair1 - pe_pair0);
    println!(
        "  engine Ewald  : {:.12}  (compute path)",
        d1.engine_ewald - d0.engine_ewald
    );
    println!(
        "  engine Ewald  : {:.12}  (energy-only path)",
        d1.energy_only - d0.energy_only
    );
    println!("  brute real    : {:.12}", d1.bf_real - d0.bf_real);
    println!("  brute recip   : {:.12}", d1.bf_rec - d0.bf_rec);
    println!(
        "  brute total   : {:.12}",
        d1.bf_total - d0.bf_total
    );
    println!();
    println!("pair PE at x0: {:.12}", pe_pair0);
    println!("pair PE at x1: {:.12}", pe_pair1);

    // --- Finite-difference force consistency check on the displaced config:
    // does the engine's Ewald force equal -grad of the engine's Ewald energy?
    // (Pinpoints whether forces and energy are mutually consistent away from
    // the symmetric lattice, where forces are nonzero.) ---
    println!();
    println!("== FD force check at x1 (engine compute() forces vs -grad energy()) ==");
    let h = 1e-6;
    let mut sysf = sys1.clone();
    sysf.zero_forces();
    let mut nlf = NeighborList::new(5.0);
    let mut vf = Virial::default();
    ewald.compute(&mut sysf, &mut vf, &mut nlf).unwrap();
    let mut worst = 0.0f64;
    let mut worst_desc = String::new();
    for i in 0..n {
        for k in 0..3 {
            let mut sp = sys1.clone();
            let p = sp.pos(i);
            let mut a = p.0;
            a[k] += h;
            sp.set_pos(i, Vec3(a));
            let mut sm = sys1.clone();
            let p = sm.pos(i);
            let mut a = p.0;
            a[k] -= h;
            sm.set_pos(i, Vec3(a));
            let ep = ewald.energy(&sp, &mut NeighborList::new(5.0)).unwrap();
            let em = ewald.energy(&sm, &mut NeighborList::new(5.0)).unwrap();
            let fd = -(ep - em) / (2.0 * h);
            let f = sysf.force(i).0[k];
            let err = (f - fd).abs();
            if err > worst {
                worst = err;
                worst_desc = format!("atom {i} axis {k}: force {f:.9e} vs FD {fd:.9e}");
            }
        }
    }
    println!("  worst |F - FD|: {worst:.3e}  ({worst_desc})");
}

// ---------------------------------------------------------------------------
// SHAKE measurement
// ---------------------------------------------------------------------------

fn shake_measure() {
    // Exactly the bench/validation/shake/modena.toml system: 16 rigid
    // water-like trimers (O + 2 H, all masses 1.0), 48 atoms, box 10^3.
    let base = Path::new("bench/validation/shake");
    let mut system =
        modena_core::read_xyz(base.join("trimers.xyz"), SimBox::from_array([10.0, 10.0, 10.0]))
            .unwrap();
    for m in system.masses.iter_mut() {
        *m = 1.0;
    }
    init_velocities(&mut system, 0.5, 87287);

    // LJ pair (sigma=1, epsilon=0.1, cutoff=2.5) + k=0 bonds/angles.
    let lj = LennardJones::new(1.0, 0.1, 2.5).unwrap();
    let pair = PairForceModel::new(lj);
    let hh = (8.0f64 / 3.0).sqrt();
    let theta = (-1.0f64 / 3.0).acos();
    let mut topo = BondedTopology::new();
    for m in 0..16 {
        let (o, h1, h2) = (3 * m, 3 * m + 1, 3 * m + 2);
        topo.bonds.push(Bond { i: o, j: h1, r0: 1.0, k: 0.0 });
        topo.bonds.push(Bond { i: o, j: h2, r0: 1.0, k: 0.0 });
        topo.angles.push(modena_core::Angle {
            i: h1,
            j: o,
            k: h2,
            theta0: theta,
            ktheta: 0.0,
        });
    }
    let bonded = BondedForceModel::new(topo);
    let model = CompositeForceModel::new(pair).with_bonded(bonded);

    // Constraints: 3 per trimer (2 O-H bonds + the H-H 1-3 distance).
    let mut cons = Constraints::new();
    for m in 0..16 {
        let (o, h1, h2) = (3 * m, 3 * m + 1, 3 * m + 2);
        cons.bonds.push(Constraint { i: o, j: h1, d: 1.0 });
        cons.bonds.push(Constraint { i: o, j: h2, d: 1.0 });
        cons.bonds.push(Constraint { i: h1, j: h2, d: hh });
    }
    let mut integrator = VelocityVerlet::with_constraints(cons.clone());

    let n = system.n();
    let dt = 0.005_f64;
    let steps = 2000usize;
    let mut e_first = 0.0f64; // E at step 1 (post initial RATTLE projection)

    // Initial force evaluation (the runner's from_config path).
    let mut neighbors = NeighborList::new(2.5);
    let mut virial = Virial::default();
    system.zero_forces();
    let pe0 = model.compute(&mut system, &mut virial, &mut neighbors).unwrap();

    let ke0 = system.kinetic_energy();
    println!("atoms: {n}");
    println!("step 0: PE = {:.12}  KE = {:.12}  E = {:.12}", pe0, ke0, pe0 + ke0);

    let bonds = &cons.bonds;
    let bond_len = |sys: &System, c: &Constraint| -> f64 {
        sys.box_size.min_image(sys.pos(c.j) - sys.pos(c.i)).norm()
    };

    let e0 = pe0 + ke0;
    let mut max_len_dev_all = 0.0f64;
    let mut sum_sq = 0.0f64;
    let mut count = 0usize;
    let mut max_rattle = 0.0f64;
    let mut e_min = f64::MAX;
    let mut e_max = f64::MIN;
    let mut e_final = 0.0f64;

    for step in 1..=steps {
        let e = integrator
            .step(&mut system, &model, &mut neighbors, dt)
            .unwrap();
        e_final = e.total;
        if step == 1 {
            e_first = e.total;
        }
        // Bond residuals at full precision.
        let mut step_max = 0.0f64;
        for c in bonds {
            let dev = (bond_len(&system, c) - c.d).abs();
            if dev > step_max {
                step_max = dev;
            }
            sum_sq += dev * dev;
            count += 1;
        }
        if step_max > max_len_dev_all {
            max_len_dev_all = step_max;
        }
        // RATTLE residual r·v.
        for c in bonds {
            let rij = system.box_size.min_image(system.pos(c.j) - system.pos(c.i));
            let vij = system.vel(c.j) - system.vel(c.i);
            let rv = rij.dot(vij).abs();
            if rv > max_rattle {
                max_rattle = rv;
            }
        }
        if step >= 100 {
            e_min = e_min.min(e.total);
            e_max = e_max.max(e.total);
        }
        if step % 200 == 0 {
            println!(
                "step {step:>4}: E = {:.10}  T = {:.8}  max|r-d| this run = {:.3e}",
                e.total,
                e.temperature,
                max_len_dev_all
            );
        }
    }
    let rms = (sum_sq / count as f64).sqrt();
    println!();
    println!("== SHAKE summary ==");
    println!("  E0 (pre-projection) : {:.12}", e0);
    println!("  E step 1 (projected): {:.12}  (one-time RATTLE projection of the initial velocities)", e_first);
    println!("  E final             : {:.12}", e_final);
    println!("  E range steps>=100  : [{:.12}, {:.12}]", e_min, e_max);
    println!("  rel. drift >=100    : {:.3e}", (e_max - e_min).abs() / e_max.abs());
    println!("  bond |r-d| RMS      : {rms:.6e}  (over {count} bond samples)");
    println!("  bond |r-d| max      : {max_len_dev_all:.6e}");
    println!("  RATTLE |r·v| max    : {max_rattle:.6e}");
}

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_default();
    match mode.as_str() {
        "ewald-decompose" => ewald_decompose(),
        "shake-measure" => shake_measure(),
        _ => {
            eprintln!("usage: validation-harness <ewald-decompose|shake-measure>");
            std::process::exit(2);
        }
    }
}
