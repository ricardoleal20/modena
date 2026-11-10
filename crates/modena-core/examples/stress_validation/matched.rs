//! The reference adapter uses a fixed roster, fixture schema and Hamiltonian.
use super::{
    ensemble::gas_fixture,
    fixtures::{analytic_endpoint, Zero},
    protocol::*,
};
use modena_core::stress::{endpoint_pressure, ConstraintSolverOptions};
use modena_core::{
    CompositeForceModel, EwaldForceModel, ForceModel, Integrator, LangevinBAOAB, LennardJones,
    LjTypeParams, MixingRule, NeighborList, PairForceModel, PairScale, SpmeForceModel, System,
    VelocityVerlet,
};
use serde::{Deserialize, Serialize};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Job {
    pub id: String,
    pub fixture: String,
    pub backend: String,
    pub dt: f64,
    pub steps: u64,
    pub burn_steps: u64,
    pub interval: u64,
    pub seed: Option<u64>,
}
pub fn roster() -> Vec<Job> {
    let mut a = vec![];
    for fixture in ["rotor", "water-nonprincipal"] {
        for n in [100, 200, 400, 800] {
            a.push(Job {
                id: format!("{fixture}-nve-{n}"),
                fixture: fixture.into(),
                backend: "zero".into(),
                dt: 8. / (0.7 * n as f64),
                steps: n,
                burn_steps: 0,
                interval: 1,
                seed: None,
            });
        }
    }
    for backend in ["ewald", "spme"] {
        for n in [200, 400, 800, 1600] {
            a.push(Job {
                id: format!("{backend}-nve-{n}"),
                fixture: "interacting32".into(),
                backend: backend.into(),
                dt: 2. / n as f64,
                steps: n,
                burn_steps: 0,
                interval: 1,
                seed: None,
            });
        }
        for seed in 2001..=2004 {
            for n in [100, 200] {
                a.push(Job {
                    id: format!("{backend}-nvt-{seed}-{n}"),
                    fixture: format!("interacting32-{seed}"),
                    backend: backend.into(),
                    dt: 1. / n as f64,
                    steps: 120 * n,
                    burn_steps: 20 * n,
                    interval: n / 10,
                    seed: Some(seed),
                });
            }
        }
    }
    a
}
pub fn fixtures() -> R<Vec<ValidationFixture>> {
    let mut f = super::fixtures::det_fixtures();
    for water in [false, true] {
        f.push(gas_fixture(water, false, 0)?);
        let mut preflight = gas_fixture(water, false, 2001)?;
        preflight.id.push_str("-2001");
        f.push(preflight);
    }
    for stream in [0, 2001, 2002, 2003, 2004] {
        f.push(gas_fixture(true, true, stream)?);
    }
    for (name, d, charges, excluded) in [
        ("pair-oo", 1.33, [0., 0.], false),
        ("pair-charge", 2.2, [-0.8, 0.4], false),
        ("pair-excluded-in", 1.6, [-0.8, 0.4], true),
        ("pair-excluded-out", 4.1, [-0.8, 0.4], true),
    ] {
        let mut s = System::new(modena_core::SimBox::cubic(12.));
        s.resize(2);
        s.set_pos(0, modena_core::Vec3::new(4., 4., 4.));
        s.set_pos(1, s.pos(0) + modena_core::Vec3::new(1., 2., 2.) * (d / 3.));
        s.charges = charges.into();
        s.masses = if name == "pair-oo" {
            vec![15.9994; 2]
        } else {
            vec![15.9994, 1.008]
        };
        let mut probe =
            ValidationFixture::from_system(name.into(), &s, modena_core::Constraints::new());
        probe.kappa = Some(0.8);
        probe.real_cutoff = Some(3.5);
        probe.reciprocal_cutoff = Some(8.);
        if excluded {
            probe.pair_scales.push((0, 1, 0.))
        }
        f.push(probe);
    }
    Ok(f)
}
pub fn model(f: &ValidationFixture, b: &str) -> R<(Box<dyn ForceModel>, NeighborList)> {
    if f.kappa.is_none() {
        if b != "zero" {
            return Err("free fixture requires zero model".into());
        }
        return Ok((Box::new(Zero), NeighborList::new(1.)));
    }
    let k = f.reciprocal_cutoff.ok_or("missing reciprocal cutoff")?;
    let scales: Vec<_> = f
        .pair_scales
        .iter()
        .map(|&(i, j, s)| PairScale {
            i,
            j,
            lj: s,
            coulomb: s,
        })
        .collect();
    let pair = PairForceModel::new(
        LennardJones::from_types(
            &[
                LjTypeParams {
                    sigma: 1.,
                    epsilon: 0.1,
                },
                LjTypeParams {
                    sigma: 1.,
                    epsilon: 0.,
                },
            ],
            MixingRule::LorentzBerthelot,
            &[],
            2.5,
        )
        .map_err(err)?,
    );
    let base = CompositeForceModel::new(pair);
    let model = match b {
        "ewald" | "ewald-10" | "ewald-12" => base.with_ewald(
            EwaldForceModel::new(
                0.8,
                3.5,
                match b {
                    "ewald-10" => 10.,
                    "ewald-12" => 12.,
                    _ => k,
                },
            )
            .map_err(err)?
            .with_pair_scales(scales.clone())
            .map_err(err)?,
        ),
        "spme" => base.with_spme(
            SpmeForceModel::new(0.8, 3.5, k)
                .map_err(err)?
                .with_pair_scales(scales.clone())
                .map_err(err)?,
        ),
        _ => return Err("unsupported matched backend".into()),
    };
    Ok((
        Box::new(model),
        NeighborList::with_pair_scales(3.5, 0.3, &scales).map_err(err)?,
    ))
}
#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Frame {
    pub step: u64,
    pub time: f64,
    pub positions: Vec<[f64; 3]>,
    pub velocities: Vec<[f64; 3]>,
}
fn frame(s: &System, step: u64, time: f64) -> R<Frame> {
    finite(&s.positions)?;
    finite(&s.velocities)?;
    finite(&[time])?;
    Ok(Frame {
        step,
        time,
        positions: (0..s.n()).map(|i| s.pos(i).0).collect(),
        velocities: (0..s.n()).map(|i| s.vel(i).0).collect(),
    })
}
pub fn evaluate(f: &ValidationFixture, b: &str, input: &Path, output: &Path) -> R<()> {
    let (model, mut nl) = model(f, b)?;
    let mut w = BufWriter::new(File::create_new(output).map_err(err)?);
    let mut count = 0;
    for l in BufReader::new(File::open(input).map_err(err)?).lines() {
        let row: Frame = serde_json::from_str(&l.map_err(err)?).map_err(err)?;
        finite(&[row.time])?;
        let mut snapshot = f.clone();
        snapshot.positions = row.positions;
        snapshot.velocities = row.velocities;
        let mut s = snapshot.system()?;
        nl.force_rebuild(&s);
        let e = model.evaluate(&mut s, &mut nl).map_err(err)?;
        finite(&[e.energy, e.scalar_virial])?;
        finite(&s.forces)?;
        let p = endpoint_pressure(&s, &e, Some(&f.constraints)).map_err(err)?;
        let analytic = if f.kappa.is_none() && s.n() <= 3 {
            Some(analytic_endpoint(&s)?)
        } else {
            None
        };
        let record = serde_json::json!({"reciprocal_cutoff":match b{"ewald-10"=>Some(10.),"ewald-12"=>Some(12.),_=>f.reciprocal_cutoff},"spme_mesh":if b=="spme"{Some(modena_core::spme::mesh_for_box(&s.box_size,f.reciprocal_cutoff.ok_or("missing cutoff")?))}else{None},"frame":frame(&s,row.step,row.time)?,"fixture":f.id,"backend":b,"potential_energy":e.energy,"forces":s.forces,"conservative_tensor":tensor(&e.tensor)?,"endpoint":tensor(&p.total)?,"analytic_force_pressure":analytic});
        serde_json::to_writer(&mut w, &record).map_err(err)?;
        writeln!(w).map_err(err)?;
        count += 1;
    }
    if count == 0 {
        return Err("empty snapshot stream".into());
    }
    w.flush().map_err(err)
}
pub fn trajectory(f: &ValidationFixture, j: &Job, output: &Path) -> R<()> {
    let mut s = f.system()?;
    let (model, mut nl) = model(f, &j.backend)?;
    let mut int: Box<dyn Integrator> = if let Some(seed) = j.seed {
        Box::new(LangevinBAOAB::with_constraints(
            1.,
            1.,
            seed,
            f.constraints.clone(),
        ))
    } else {
        Box::new(VelocityVerlet::with_constraints(f.constraints.clone()))
    };
    let mut prev = model.evaluate(&mut s, &mut nl).map_err(err)?;
    let mut w = BufWriter::new(File::create_new(output).map_err(err)?);
    let p0 = endpoint_pressure(&s, &prev, Some(&f.constraints)).map_err(err)?;
    finite(&[prev.energy])?;
    finite(&s.forces)?;
    let row = serde_json::json!({"frame":frame(&s,0,0.)?,"job":j.id,"potential_energy":prev.energy,"forces":s.forces,"conservative_tensor":tensor(&prev.tensor)?,"endpoint":tensor(&p0.total)?,"interval":null});
    finite(&[prev.energy])?;
    finite(&s.forces)?;
    serde_json::to_writer(&mut w, &row).map_err(err)?;
    writeln!(w).map_err(err)?;
    let mut work = Work::begin();
    let molecules = if s.n() == 96 { 32 } else { 1 };
    let mut exposure = 0;
    let mut window = [0.; 6];
    for step in 1..=j.steps {
        let t = int
            .step_with_stress_options(
                &mut s,
                &*model,
                &mut nl,
                j.dt,
                &prev,
                ConstraintSolverOptions::TIGHT,
            )
            .map_err(err)?
            .stress
            .ok_or("missing step stress")?;
        work.step(molecules)?;
        prev = t.final_model;
        let interval = tensor(&t.interval.total)?;
        if step <= j.burn_steps {
            continue;
        }
        exposure += 1;
        for k in 0..6 {
            window[k] += interval[k]
        }
        if !(step - j.burn_steps).is_multiple_of(j.interval) {
            continue;
        }
        finite(&window)?;
        finite(&s.forces)?;
        finite(&[prev.energy])?;
        let sample = if s.n() == 96 {
            Some(super::ensemble::sample(
                &s,
                3,
                j.seed.unwrap_or(0),
                step,
                j.dt,
                tensor(&t.endpoint.total)?,
                interval,
            )?)
        } else {
            None
        };
        let row = serde_json::json!({"frame":frame(&s,step,step as f64*j.dt)?,"job":j.id,"potential_energy":prev.energy,"forces":s.forces,"conservative_tensor":tensor(&prev.tensor)?,"endpoint":tensor(&t.endpoint.total)?,"interval":interval,"window_steps":exposure,"window_interval":window.map(|x|x/exposure as f64),"sample":sample});
        serde_json::to_writer(&mut w, &row).map_err(err)?;
        writeln!(w).map_err(err)?;
        w.flush().map_err(err)?;
        window = [0.; 6];
        exposure = 0;
    }
    work.persist(true)
}
