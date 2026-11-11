//! One opt-in, epoch-bound real Cu4 EMT measurement. No provider runs at preparation.
use super::protocol::{err, finite, read, tensor, R};
use modena_core::units::UnitFactors;
use modena_core::{AseEmtForceModel, AseEmtSpec, ForceModel, NeighborList, SimBox, System, Vec3};
use serde::Deserialize;
use std::{fs::File, io::Write, path::Path};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Snapshot {
    component: Option<usize>,
    eta: f64,
    cell: [[f64; 3]; 3],
    positions: [[f64; 3]; 4],
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Inputs {
    provider: AseEmtSpec,
    mass_amu: f64,
    units: [f64; 3],
    snapshots: Vec<Snapshot>,
}
pub fn measure(dir: &Path) -> R<()> {
    let input: Inputs = read(&dir.join("reference/inputs/ase.json"))?;
    if std::env::var("MODENA_ASE_PYTHON").map_err(err)? != input.provider.python
        || input.snapshots.len() != 49
        || input.provider.atomic_numbers != [29; 4]
    {
        return Err("ASE gate requires the pinned Cu4 protocol".into());
    }
    let mut output = File::create_new(dir.join("reference/ase.jsonl")).map_err(err)?;
    let mut requests = File::create_new(dir.join("reference/ase-requests.jsonl")).map_err(err)?;
    let mut request = |id: usize, kind: &str, stage: &str| -> R<()> {
        writeln!(
            requests,
            "{}",
            serde_json::json!({"id":id,"kind":kind,"stage":stage})
        )
        .map_err(err)?;
        requests.flush().map_err(err)
    };
    request(1, "info", "attempt")?;
    let [sigma, epsilon, mass] = input.units;
    let model = AseEmtForceModel::new(input.provider, UnitFactors::new(sigma, epsilon, mass))
        .map_err(err)?;
    request(1, "info", "complete")?;
    writeln!(
        output,
        "{}",
        serde_json::json!({"completed_info_requests":1,"identity":model.identity()})
    )
    .map_err(err)?;
    for (i, snapshot) in input.snapshots.into_iter().enumerate() {
        let mut system = System::new(SimBox::triclinic(snapshot.cell.map(Vec3)));
        system.resize(4);
        system.positions = snapshot.positions.into_iter().flatten().collect();
        system.masses.fill(input.mass_amu / mass);
        system.atomic_numbers.fill(29);
        // One direct wire response at the base is independent of the second request
        // through the trait's reduced-unit tensor conversion.
        let wire = if i == 0 {
            request(2, "wire-base", "attempt")?;
            let wire = model.evaluate(&system).map_err(err)?;
            request(2, "wire-base", "complete")?;
            Some(wire)
        } else {
            None
        };
        request(i + 3, "reduced", "attempt")?;
        let e =
            ForceModel::evaluate(&model, &mut system, &mut NeighborList::new(1.)).map_err(err)?;
        request(i + 3, "reduced", "complete")?;
        finite(&[e.energy, e.scalar_virial])?;
        finite(&system.forces)?;
        let row = serde_json::json!({"index":i,"completed_evaluate_requests":i+2,
            "component":snapshot.component,"eta":snapshot.eta,"cell":snapshot.cell,
            "positions":snapshot.positions,"wire":wire,"energy":e.energy,"forces":system.forces,
            "virial":tensor(&e.tensor)?,"trace":e.scalar_virial,"volume":system.box_size.volume()});
        writeln!(output, "{row}").map_err(err)?;
        output.flush().map_err(err)?;
    }
    Ok(())
}
