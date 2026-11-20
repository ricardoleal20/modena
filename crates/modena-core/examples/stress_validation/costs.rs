use super::{fixtures::Zero, protocol::*};
use modena_core::{
    constraint_stress::{constraint_geometry, endpoint_reaction, ConstraintImpulse},
    io::restart::{self, Checkpoint, ContinuationState, RestartStress, StressCheckpoint},
    stress::{ConstraintSolverOptions, ModelEvaluation},
    ForceModel, Integrator, NeighborList, PairScale, SpmeForceModel, VelocityVerlet,
};
use serde_json::json;
use std::{path::Path, time::Instant};
pub fn measure(out: &Path, p: &ValidationProtocol) -> R<()> {
    let mut work = Work::begin();
    let mut rows = vec![];
    for (id, steps) in [("rotor", 128), ("gas-water", 32)] {
        let f = p
            .fixtures
            .iter()
            .find(|f| f.id == id)
            .ok_or("cost fixture")?;
        let mut s = f.system()?;
        let mut vv = VelocityVerlet::with_constraints(f.constraints.clone());
        let mut nl = NeighborList::new(1.);
        let mut actual = 0;
        let mut last = None;
        let dt = 0.005;
        let start = Instant::now();
        for _ in 0..steps {
            last = vv
                .step_with_stress(&mut s, &Zero, &mut nl, dt, &ModelEvaluation::zero())
                .map_err(err)?
                .stress;
            work.step(if id == "rotor" { 1 } else { 32 })?;
            actual += 1
        }
        let ns = start.elapsed().as_nanos();
        let st = last.ok_or("missing cost step")?;
        let geometry = constraint_geometry(&s, &f.constraints).map_err(err)?;
        let stamp = Instant::now();
        let reaction = endpoint_reaction(&s, &geometry).map_err(err)?;
        let solve_ns = stamp.elapsed().as_nanos();
        let stage:Vec<_>=st.stages.iter().map(|x|json!({"stage":x.stage,"entries":x.impulses.len(),"live_entry_bytes":x.impulses.len()*std::mem::size_of::<ConstraintImpulse>(),"capacity_bytes":x.impulses.capacity()*std::mem::size_of::<ConstraintImpulse>(),"iterations":x.iterations})).collect();
        let sig =
            json!({"integrator":{"type":"velocity_verlet","dt":dt},"constraints":f.constraints});
        let cp = Checkpoint {
            atoms: s.n(),
            step: actual,
            time: actual as f64 * dt,
            box_size: s.box_size,
            positions: s.positions.clone(),
            velocities: s.velocities.clone(),
            continuation: Some(ContinuationState {
                forces: s.forces.clone(),
                masses: s.masses.clone(),
                species: s.species.clone(),
                atomic_numbers: s.atomic_numbers.clone(),
                charges: s.charges.clone(),
                elements: s.elements.clone(),
                constraint_count: f.constraints.count(),
                conserved_com_dof: 0,
                potential_energy: 0.,
                virial: 0.,
                physical_signature: sig.clone(),
                integrator: vv.checkpoint_state().map_err(err)?,
                thermostat: serde_json::Value::Null,
                barostat: serde_json::Value::Null,
            }),
        };
        let mut signature = sig;
        signature["pressure_estimator"] = json!({"version":1,"solver":ConstraintSolverOptions::DEFAULT,"endpoint_timing":"endpoint_after_controllers","interval_timing":"interval_pre_controllers"});
        let v2 = out.join(format!("cost-{id}-v2.chk"));
        let v3 = out.join(format!("cost-{id}-v3.chk"));
        restart::save(&v2, &cp).map_err(err)?;
        restart::save_with_stress(
            &v3,
            &StressCheckpoint {
                mechanical: cp,
                stress: Some(RestartStress {
                    estimator_version: 1,
                    solver: st.solver,
                    model: st.final_model.clone(),
                    endpoint: st.endpoint.clone(),
                    interval: Some(st.interval.clone()),
                    last_stages: st.stages.clone(),
                    last_bath: None,
                    step: actual,
                    dt,
                    physical_signature: signature,
                }),
            },
        )
        .map_err(err)?;
        let pressure_path = out.join(format!("cost-{id}-pressure.jsonl"));
        #[derive(serde::Serialize)]
        struct OutputRow<'a> {
            step: usize,
            time: f64,
            endpoint: &'a modena_core::stress::PressureObservation,
            interval: Option<&'a modena_core::stress::IntervalPressure>,
        }
        write_new(
            &pressure_path,
            &OutputRow {
                step: actual,
                time: actual as f64 * dt,
                endpoint: &st.endpoint,
                interval: Some(&st.interval),
            },
        )?;
        let csv = out.join(format!("cost-{id}-energy.csv"));
        let mut log = modena_core::LogWriter::create(&csv).map_err(err)?;
        let e = modena_core::Energies::from(&s, 0., &modena_core::Virial::default());
        log.log(
            actual,
            actual as f64 * dt,
            e.kinetic,
            e.potential,
            e.total,
            e.temperature,
            e.pressure,
            e.virial,
        )
        .map_err(err)?;
        log.flush().map_err(err)?;
        let mut bath_state = f.system()?;
        let mut langevin =
            modena_core::LangevinBAOAB::with_constraints(1., 1., 1001, f.constraints.clone());
        let bath = langevin
            .step_with_stress_options(
                &mut bath_state,
                &Zero,
                &mut nl,
                dt,
                &ModelEvaluation::zero(),
                ConstraintSolverOptions::TIGHT,
            )
            .map_err(err)?
            .stress
            .ok_or("missing BAOAB cost stage")?;
        work.step(if id == "rotor" { 1 } else { 32 })?;
        let bath_stages:Vec<_>=bath.stages.iter().map(|r|json!({"stage":r.stage,"entries":r.impulses.len(),"live_entry_bytes":r.impulses.len()*std::mem::size_of::<ConstraintImpulse>(),"capacity_bytes":r.impulses.capacity()*std::mem::size_of::<ConstraintImpulse>()})).collect();
        rows.push(json!({"fixture":id,"pressure_output_jsonl_bytes":std::fs::metadata(pressure_path).map_err(err)?.len(),"energy_csv_including_header_bytes":std::fs::metadata(csv).map_err(err)?.len(),"baoab_stages":bath_stages,"baoab_bath_json_bytes":serde_json::to_vec(&bath.bath).map_err(err)?.len(),"reporting_whole_step_ns":ns,"actual_steps":actual,"dt":dt,"final_q":s.positions,"final_v":s.velocities,"stages":stage,"rollback_bytes_per_projection":24*s.n(),"reconstruction_bytes_per_reporting_projection":48*s.n(),"block_edge_counts":if id=="rotor"{vec![1]}else{vec![3;32]},"endpoint_reaction_ns_including_validation":solve_ns,"relative_linear_residual":reaction.relative_linear_residual,"step_stress_json_bytes":serde_json::to_vec(&st).map_err(err)?.len(),"v2_bytes":std::fs::metadata(v2).map_err(err)?.len(),"v3_bytes":std::fs::metadata(v3).map_err(err)?.len()}));
    }
    let f = p
        .fixtures
        .iter()
        .find(|f| f.id == "interacting32")
        .ok_or("profile fixture")?;
    let mut s = f.system()?;
    let scales: Vec<_> = f
        .pair_scales
        .iter()
        .map(|&(i, j, v)| PairScale {
            i,
            j,
            lj: v,
            coulomb: v,
        })
        .collect();
    let model = SpmeForceModel::new(0.8, 3.5, 8.)
        .map_err(err)?
        .with_pair_scales(scales.clone())
        .map_err(err)?;
    let mut nl = NeighborList::with_pair_scales(3.5, 0.3, &scales).map_err(err)?;
    let mesh = modena_core::spme::mesh_for_box(&s.box_size, 8.);
    let mut profile = vec![];
    model.evaluate(&mut s, &mut nl).map_err(err)?;
    for pair in 0..8 {
        for mode in [pair % 2 == 0, pair % 2 != 0] {
            s.forces.fill(0.);
            let start = Instant::now();
            let (ev, p) = if mode {
                let (v, p) = model.evaluate_profiled(&mut s, &mut nl).map_err(err)?;
                (v, Some(p.total.as_nanos()))
            } else {
                (model.evaluate(&mut s, &mut nl).map_err(err)?, None)
            };
            let elapsed = start.elapsed().as_nanos();
            finite(&s.forces)?;
            finite(&[ev.energy])?;
            profile.push(json!({"profiled":mode,"elapsed_ns":elapsed,"inside_profile_ns":p,"energy":ev.energy,"tensor":tensor(&ev.tensor)?}));
        }
    }
    work.persist(true)?;
    write_new(
        &out.join("costs.json"),
        &json!({"steps":rows,"observer_pairs":profile,"profile_warmup_calls":1,"mesh":mesh,"workspace_factor_delta_bytes":48*mesh.iter().product::<usize>()}),
    )
}
