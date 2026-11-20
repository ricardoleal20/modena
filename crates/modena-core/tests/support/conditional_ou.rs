//! Unit tests of the same fixed-q OU operator used between both BAOAB drifts.
use super::projected_ou;
use crate::constraint_stress::ProjectionStage;
use crate::stress::ConstraintSolverOptions;
use crate::{Constraints, System, Vec3};
use chacha20::ChaCha12Rng;
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};

#[allow(dead_code)]
mod fixtures {
    use crate as modena_core;
    include!("stress.rs");
}

const M: [f64; 3] = [16., 1., 2.];
const DIM: usize = 9;
type Vector = [f64; DIM];
type Matrix = [Vector; DIM];

fn water() -> (System, Constraints) {
    let (mut system, cons) = fixtures::water(M);
    let omega = Vec3::new(0.2, 0.3, 0.4);
    for i in 0..3 {
        let rho = system.pos(i) - Vec3::new(5., 6., 7.);
        system.set_vel(i, Vec3::new(0.2, -0.1, 0.3) + omega.cross(rho));
    }
    (system, cons)
}

#[test]
fn fixture_uses_global_nonprincipal_angular_velocity() {
    let (system, _) = water();
    // The coordinates already carry the frozen Rodrigues rotation. The angular
    // velocity in the specification is expressed in the laboratory frame.
    for i in 0..system.n() {
        let [x, y, z] = (system.pos(i) - Vec3::new(5., 6., 7.)).0;
        let expected = [
            0.2 + 0.3 * z - 0.4 * y,
            -0.1 + 0.4 * x - 0.2 * z,
            0.3 + 0.2 * y - 0.3 * x,
        ];
        for (actual, target) in system.vel(i).0.into_iter().zip(expected) {
            fixtures::near(actual, target, 1e-14);
        }
    }
}

/// Form the full orthogonal projector on the kernel of G M^{-1/2}.
/// The three rows come directly from the triangle coordinates. Orthonormalizing
/// them avoids sharing any production Gram solver, edge indexing, or ledger.
fn covariance(system: &System, c: f64) -> Matrix {
    let mut basis = [[0.; DIM]; 3];
    for (k, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
        let r = system.pos(j) - system.pos(i);
        let mut row = [0.; DIM];
        for axis in 0..3 {
            row[3 * i + axis] = -r.0[axis] / M[i].sqrt();
            row[3 * j + axis] = r.0[axis] / M[j].sqrt();
        }
        for unit in &basis[..k] {
            let dot: f64 = unit.iter().zip(row).map(|(u, v)| u * v).sum();
            for (x, u) in row.iter_mut().zip(unit) {
                *x -= dot * u;
            }
        }
        let norm = row.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(norm > 1e-10);
        basis[k] = row.map(|x| x / norm);
    }
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let normal: f64 = basis.iter().map(|row| row[i] * row[j]).sum();
            (1. - c * c) * (M[i / 3] * M[j / 3]).sqrt() * (f64::from(i == j) - normal)
        })
    })
}

#[derive(Default)]
struct Moments {
    count: usize,
    sums: Vector,
    products: Matrix,
}
impl Moments {
    fn observe(&mut self, x: Vector) {
        self.count += 1;
        for (i, &xi) in x.iter().enumerate() {
            self.sums[i] += xi;
            for (j, &xj) in x.iter().enumerate() {
                self.products[i][j] += xi * xj;
            }
        }
    }
    fn mean(&self) -> Vector {
        self.sums.map(|s| s / self.count as f64)
    }
    fn covariance(&self) -> Matrix {
        let n = self.count as f64;
        std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (self.products[i][j] - self.sums[i] * self.sums[j] / n) / (n - 1.)
            })
        })
    }
}

#[test]
fn fixed_water_conditional_momentum_law() {
    let (mut system, cons) = water();
    let initial = system.velocities.clone();
    let c = (-0.2_f64).exp();
    let expected = covariance(&system, c);
    let expected_raw: Vector = std::array::from_fn(|i| (c - 1.) * initial[i] * M[i / 3]);
    let mut raw = Moments::default();
    let mut projected = Moments::default();
    let mut reaction = Moments::default();
    let mut counts = [0usize; 8];
    for (seed, count) in counts.iter_mut().enumerate() {
        let mut rng = ChaCha12Rng::seed_from_u64(seed as u64);
        for _ in 0..8192 {
            system.velocities.copy_from_slice(&initial);
            let (bath, report) = projected_ou(
                &mut system,
                Some(&cons),
                &mut rng,
                c,
                1.,
                ConstraintSolverOptions::TIGHT,
                true,
            )
            .unwrap()
            .unwrap();
            assert_eq!(report.stage, ProjectionStage::BaoabO);
            // Measure outgoing momentum directly from the physical state.
            projected.observe(std::array::from_fn(|i| system.velocities[i] * M[i / 3]));
            raw.observe(std::array::from_fn(|i| bath.unprojected[i / 3].0[i % 3]));
            reaction.observe(std::array::from_fn(|i| bath.reaction[i / 3].0[i % 3]));
            *count += 1;
        }
        println!("conditional seed={seed} draws={count}");
    }
    let n = projected.count as f64;
    assert_eq!(counts, [8192; 8]);
    assert_eq!(projected.count, 65536);
    let cov = projected.covariance();
    let pm = projected.mean();
    let rm = reaction.mean();
    let um = raw.mean();
    let expected_p: Vector = std::array::from_fn(|i| c * initial[i] * M[i / 3]);
    let normalized: Matrix = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (cov[i][j] - expected[i][j]) / (expected[i][i] * expected[j][j]).sqrt()
        })
    });
    let se: Matrix = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            ((expected[i][i] * expected[j][j] + expected[i][j].powi(2)) / (n - 1.)).sqrt()
        })
    });
    println!(
        "conditional_statistics {}",
        serde_json::json!({
            "positions":system.positions, "masses":M, "velocities":initial, "decay":c, "kT":1,
            "draws_by_seed":counts, "total":projected.count, "raw_mean":um,"expected_raw_mean":expected_raw,
            "reaction_mean":rm,"projected_mean":pm,"expected_projected_mean":expected_p,
            "raw_covariance":raw.covariance(),"reaction_covariance":reaction.covariance(),
            "projected_covariance":cov,"expected_projected_covariance":expected,
            "normalized_covariance_error":normalized,"covariance_se":se
        })
    );
    for i in 0..DIM {
        let raw_variance = (1. - c * c) * M[i / 3];
        let reaction_variance = raw_variance - expected[i][i];
        for (name, error, variance) in [
            ("projected", pm[i] - expected_p[i], expected[i][i]),
            ("raw", um[i] - expected_raw[i], raw_variance),
            ("reaction", rm[i], reaction_variance),
        ] {
            let limit = 4. * (variance / n).sqrt() + 1e-3;
            assert!(
                error.abs() <= limit,
                "{name} mean [{i}] error={error}, limit={limit}"
            );
        }
        for j in 0..DIM {
            let error = (cov[i][j] - expected[i][j]).abs();
            assert!(
                normalized[i][j].abs() <= 0.02 && error <= 4. * se[i][j],
                "covariance [{i},{j}] error={error}, normalized={}, 4SE={}",
                normalized[i][j],
                4. * se[i][j]
            );
        }
    }
}

#[test]
fn zero_amplitude_limits_consume_all_normal_draws() {
    for (gamma, h) in [(0., 0.2_f64), (1., 0.)] {
        let (mut s, cons) = water();
        let initial = s.velocities.clone();
        let mut rng = ChaCha12Rng::seed_from_u64(1001);
        let mut reference = ChaCha12Rng::seed_from_u64(1001);
        let (bath, projection) = projected_ou(
            &mut s,
            Some(&cons),
            &mut rng,
            (-gamma * h).exp(),
            1.,
            ConstraintSolverOptions::TIGHT,
            true,
        )
        .unwrap()
        .unwrap();
        assert_eq!(s.velocities, initial);
        for impulses in [bath.unprojected, bath.reaction] {
            assert_eq!(impulses, vec![Vec3::zero(); 3]);
        }
        assert!(projection.impulses.is_empty());
        for _ in 0..DIM {
            let _: f64 = StandardNormal.sample(&mut reference);
        }
        assert_eq!(rng.get_word_pos(), reference.get_word_pos());
        assert_eq!(rng.get_stream(), reference.get_stream());
    }
}
