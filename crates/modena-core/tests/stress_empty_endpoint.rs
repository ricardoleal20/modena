//! Isolated single-test executable so the allocator measures only this call.
#[allow(dead_code)]
#[path = "support/pair_scratch_alloc.rs"]
mod alloc;

use modena_core::constraint_stress::{constraint_geometry, endpoint_reaction};
use modena_core::stress::{endpoint_pressure, ModelEvaluation, SymTensor6};
use modena_core::{Constraints, SimBox, System, Vec3};
use std::sync::atomic::Ordering;

#[global_allocator]
static ALLOCATOR: alloc::CountingAllocator = alloc::CountingAllocator;

#[test]
fn empty_pressure_checks_state_without_allocating_constraint_work() {
    let mut s = System::new(SimBox::cubic(10.));
    s.resize(64);
    let model = ModelEvaluation::zero();
    let empty = Constraints::new();
    let geometry = constraint_geometry(&s, &empty).unwrap();
    let full = endpoint_reaction(&s, &geometry).unwrap();
    assert_eq!(full.forces, vec![Vec3::zero(); 64]);
    assert!(full.multipliers.is_empty());
    for constraints in [None, Some(&empty)] {
        alloc::reset(0);
        let live = alloc::LIVE_BYTES.load(Ordering::SeqCst);
        let p = endpoint_pressure(&s, &model, constraints).unwrap();
        let extra = alloc::HIGH_WATER.load(Ordering::SeqCst) - live;
        assert_eq!(p.constraint, SymTensor6::ZERO);
        assert_eq!(
            extra, 0,
            "empty pressure allocated {extra} additional live bytes"
        );
        for field in ["positions", "velocities", "forces", "shape", "mass", "cell"] {
            let mut bad = s.clone();
            match field {
                "positions" => bad.positions[1] = f64::NAN,
                "velocities" => bad.velocities[1] = f64::INFINITY,
                "forces" => bad.forces[1] = f64::NAN,
                "shape" => {
                    bad.velocities.pop();
                }
                "mass" => bad.masses[1] = -1.,
                _ => bad.box_size = SimBox::from_array([-10., -10., 10.]),
            }
            assert!(
                endpoint_pressure(&bad, &model, constraints).is_err(),
                "{field}"
            );
        }
    }
}
