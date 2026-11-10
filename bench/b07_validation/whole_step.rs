//! This exact driver is built against each source archive and runs the common legacy VV API.
use modena_core::{
    Constraints, ForceModel, Integrator, NeighborList, SimBox, System, VelocityVerlet, Virial,
};
use std::{fs::File, path::Path, time::Instant};
struct Zero;
impl ForceModel for Zero {
    fn compute(
        &self,
        _: &mut System,
        _: &mut Virial,
        _: &mut NeighborList,
    ) -> modena_core::Result<f64> {
        Ok(0.)
    }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let p = std::env::args().nth(1).ok_or("fixture path")?;
    let f: serde_json::Value = serde_json::from_reader(File::open(Path::new(&p))?)?;
    let steps = match f["id"].as_str() {
        Some("rotor") => 128,
        Some("gas-water") => 32,
        _ => return Err("unregistered benchmark".into()),
    };
    let positions: Vec<[f64; 3]> = serde_json::from_value(f["positions"].clone())?;
    let velocities: Vec<[f64; 3]> = serde_json::from_value(f["velocities"].clone())?;
    let cell: [[f64; 3]; 3] = serde_json::from_value(f["cell"].clone())?;
    let mut s = System::new(SimBox::from_array([cell[0][0], cell[1][1], cell[2][2]]));
    s.resize(positions.len());
    s.positions = positions.into_iter().flatten().collect();
    s.velocities = velocities.into_iter().flatten().collect();
    s.masses = serde_json::from_value(f["masses"].clone())?;
    let c: Constraints = serde_json::from_value(f["constraints"].clone())?;
    let mut vv = VelocityVerlet::with_constraints(c);
    let mut nl = NeighborList::new(1.);
    let mut actual = 0;
    let start = Instant::now();
    for _ in 0..steps {
        std::hint::black_box(vv.step(&mut s, &Zero, &mut nl, 0.005)?);
        actual += 1;
    }
    let ns = start.elapsed().as_nanos();
    if s.positions
        .iter()
        .chain(&s.velocities)
        .any(|x| !x.is_finite())
    {
        return Err("nonfinite cost state".into());
    }
    println!(
        "{}",
        serde_json::json!({"fixture":f["id"],"steps":actual,"dt":0.005,"integrator":"VV/common default constraints","elapsed_ns":ns,"final_q":s.positions,"final_v":s.velocities})
    );
    Ok(())
}
