//! Rigid body fixtures and independent laboratory-frame force / body-frame motion oracles.
use super::protocol::*;
use modena_core::constraint_stress::{constraint_geometry, endpoint_reaction};
use modena_core::stress::{ConstraintSolverOptions, ModelEvaluation};
use modena_core::{
    Constraints, ForceModel, Integrator, LangevinBAOAB, NeighborList, System, Vec3, VelocityVerlet,
    Virial,
};
use std::path::Path;
#[allow(dead_code)]
#[path = "../../tests/support/stress.rs"]
pub mod common;

pub const CENTER: Vec3 = Vec3::new(5., 6., 7.);
pub const V: Vec3 = Vec3::new(0.2, -0.1, 0.3);
pub const COMPONENTS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
pub struct Zero;
impl ForceModel for Zero {
    fn compute(
        &self,
        _: &mut System,
        _: &mut Virial,
        _: &mut NeighborList,
    ) -> modena_core::Result<f64> {
        Ok(0.)
    }
    fn evaluate(
        &self,
        _: &mut System,
        _: &mut NeighborList,
    ) -> modena_core::Result<ModelEvaluation> {
        Ok(ModelEvaluation::zero())
    }
}
pub fn det_fixtures() -> Vec<ValidationFixture> {
    let mut out = vec![];
    for (id, m1, m2) in [("rotor", 1., 3.), ("rotor-heavy", 2., 5.)] {
        let (s, c) = common::rotor(m1, m2);
        out.push(ValidationFixture::from_system(id.into(), &s, c));
    }
    for (suffix, m) in [("", [15.9994, 1.008, 1.008]), ("-heavy", [16., 1., 2.])] {
        for kind in ["principal", "nonprincipal"] {
            let (mut s, c) = common::water(m);
            if kind == "nonprincipal" {
                for i in 0..3 {
                    s.set_vel(i, V + Vec3::new(0.2, 0.3, 0.4).cross(s.pos(i) - CENTER))
                }
            }
            out.push(ValidationFixture::from_system(
                format!("water-{kind}{suffix}"),
                &s,
                c,
            ));
        }
    }
    out
}
pub fn dyad(a: Vec3, b: Vec3) -> [f64; 6] {
    COMPONENTS.map(|(i, j)| 0.5 * (a.0[i] * b.0[j] + a.0[j] * b.0[i]))
}
fn com(x: &[Vec3], m: &[f64]) -> Vec3 {
    x.iter().zip(m).fold(Vec3::zero(), |a, (&v, m)| a + v * *m) * (1. / m.iter().sum::<f64>())
}
pub fn inertia(r: &[Vec3], m: &[f64], v: Vec3) -> Vec3 {
    r.iter()
        .zip(m)
        .fold(Vec3::zero(), |a, (&x, m)| a + x.cross(v.cross(x)) * *m)
}
pub fn invert(columns: [Vec3; 3], rhs: Vec3) -> R<Vec3> {
    let [a, b, c] = columns;
    let d = a.dot(b.cross(c));
    if !d.is_finite() || d <= 0. {
        return Err("singular Cartesian inertia".into());
    }
    let v = Vec3::new(
        rhs.dot(b.cross(c)),
        rhs.dot(c.cross(a)),
        rhs.dot(a.cross(b)),
    ) * (1. / d);
    finite(&v.0)?;
    Ok(v)
}
fn axes() -> [Vec3; 3] {
    [
        Vec3::new(1., 0., 0.),
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 0., 1.),
    ]
}
pub fn relative(s: &System, start: usize, n: usize) -> (Vec<Vec3>, Vec<Vec3>, Vec3) {
    let p: Vec<_> = (start..start + n)
        .map(|i| s.box_size.min_image(s.pos(i) - s.pos(start)))
        .collect();
    let v: Vec<_> = (start..start + n).map(|i| s.vel(i)).collect();
    let m = &s.masses[start..start + n];
    let (x, u) = (com(&p, m), com(&v, m));
    (
        p.iter().map(|&q| q - x).collect(),
        v.iter().map(|&v| v - u).collect(),
        u,
    )
}
/// Compute accelerations from angular momentum in Cartesian space, independently of edge solves.
pub fn analytic_endpoint(s: &System) -> R<(Vec<Vec3>, [f64; 6])> {
    let (r, v, _) = relative(s, 0, s.n());
    let m = &s.masses;
    let f: Vec<Vec3> = match s.n() {
        2 => {
            let d = r[1] - r[0];
            let dv = v[1] - v[0];
            let mu = m[0] * m[1] / (m[0] + m[1]);
            let lam = mu
                * (dv.norm2() + d.dot(s.force(1) * (1. / m[1]) - s.force(0) * (1. / m[0])))
                / d.norm2();
            vec![d * lam, -d * lam]
        }
        3 => {
            let columns = axes().map(|a| inertia(&r, m, a));
            let l = (0..3).fold(Vec3::zero(), |a, i| a + r[i].cross(v[i]) * m[i]);
            let w = invert(columns, l)?;
            let dw = invert(columns, -w.cross(l))?;
            (0..3)
                .map(|i| (dw.cross(r[i]) + w.cross(w.cross(r[i]))) * m[i])
                .collect()
        }
        _ => return Err("analytic endpoint requires a single free rigid molecule".into()),
    };
    for v in &f {
        finite(&v.0)?
    }
    let p = std::array::from_fn(|k| {
        (0..s.n())
            .map(|i| dyad(s.vel(i), s.vel(i))[k] * m[i] + dyad(r[i], f[i])[k])
            .sum::<f64>()
            / s.box_size.volume()
    });
    finite(&p)?;
    Ok((f, p))
}
#[derive(Copy, Clone)]
struct Motion([f64; 7]);
impl Motion {
    fn add(self, b: Self, h: f64) -> Self {
        Self(std::array::from_fn(|i| self.0[i] + h * b.0[i]))
    }
    fn derivative(self, cols: [Vec3; 3]) -> R<Self> {
        let [q, x, y, z, a, b, c] = self.0;
        let w = Vec3::new(a, b, c);
        let l = cols[0] * a + cols[1] * b + cols[2] * c;
        let alpha = invert(cols, -w.cross(l))?;
        Ok(Self([
            (-x * a - y * b - z * c) * 0.5,
            (q * a + y * c - z * b) * 0.5,
            (q * b + z * a - x * c) * 0.5,
            (q * c + x * b - y * a) * 0.5,
            alpha.0[0],
            alpha.0[1],
            alpha.0[2],
        ]))
    }
}
pub fn rotate(q: [f64; 4], v: Vec3) -> Vec3 {
    let u = Vec3::new(q[1], q[2], q[3]);
    v + u.cross(v) * (2. * q[0]) + u.cross(u.cross(v)) * 2.
}
fn reference(s: &System, id: &str, n: usize) -> R<System> {
    let time = 8. / 0.7;
    let r: Vec<_> = (0..s.n()).map(|i| s.pos(i) - CENTER).collect();
    let w = if id.starts_with("rotor") {
        (s.pos(1) - s.pos(0)).cross(s.vel(1) - s.vel(0)) * (1. / 1.3_f64.powi(2))
    } else if id.contains("nonprincipal") {
        Vec3::new(0.2, 0.3, 0.4)
    } else {
        common::water_rotation(Vec3::new(0., 0., 0.7))
    };
    let (q, angular) = if id.contains("nonprincipal") {
        let cols = axes().map(|a| inertia(&r, &s.masses, a));
        let h = time / n as f64;
        let mut z = Motion([1., 0., 0., 0., w.0[0], w.0[1], w.0[2]]);
        for _ in 0..n {
            let k1 = z.derivative(cols)?;
            let k2 = z.add(k1, h / 2.).derivative(cols)?;
            let k3 = z.add(k2, h / 2.).derivative(cols)?;
            let k4 = z.add(k3, h).derivative(cols)?;
            z = z
                .add(k1, h / 6.)
                .add(k2, h / 3.)
                .add(k3, h / 3.)
                .add(k4, h / 6.);
        }
        let [a, b, c, d, x, y, z] = z.0;
        ([a, b, c, d], Vec3::new(x, y, z))
    } else {
        let axis = w * (1. / w.norm());
        let half = w.norm() * time / 2.;
        let u = axis * half.sin();
        ([half.cos(), u.0[0], u.0[1], u.0[2]], w)
    };
    let mut end = s.clone();
    for (i, &r) in r.iter().enumerate() {
        end.set_pos(i, CENTER + V * time + rotate(q, r));
        end.set_vel(i, V + rotate(q, angular.cross(r)));
    }
    finite(&end.positions)?;
    finite(&end.velocities)?;
    Ok(end)
}
fn moment(s: &System) -> [f64; 6] {
    std::array::from_fn(|k| {
        (0..s.n())
            .map(|i| dyad(s.pos(i), s.vel(i))[k] * s.masses[i])
            .sum()
    })
}
fn maxdiff(a: &[f64], b: &[f64]) -> R<f64> {
    finite(a)?;
    finite(b)?;
    if a.len() != b.len() || a.is_empty() {
        return Err("inconsistent comparison shape".into());
    }
    let d: Vec<_> = a.iter().zip(b).map(|(x, y)| (x - y).abs()).collect();
    finite(&d)?;
    Ok(d.into_iter().fold(0., f64::max))
}
fn state_err(s: &System, exact: &System) -> R<[f64; 2]> {
    Ok([
        maxdiff(&s.positions, &exact.positions)?,
        maxdiff(&s.velocities, &exact.velocities)?,
    ])
}
fn finite_maps(
    g: &mut Vec<GateResult>,
    s0: &System,
    c: &Constraints,
    id: &str,
    work: &mut Work,
    raw: &mut impl std::io::Write,
) -> R<()> {
    let h = 0.08 / 0.7;
    let m = &s0.masses;
    let mass = m[0] + m[1];
    let mu = m[0] * m[1] / mass;
    for half in [false, true] {
        let n = if half { 2 } else { 1 };
        let tau = h / n as f64;
        let a = 0.7 * tau;
        let cs = (1. - a * a).sqrt();
        let mut errors = vec![0_f64; n * 2 + 3];
        for x in [-1, 0, 1] {
            for y in [-1, 0, 1] {
                for z in [-1, 0, 1] {
                    for shift in [Vec3::zero(), Vec3::new(0.31, -0.27, 0.19)] {
                        let center = Vec3(std::array::from_fn(|i| match [x, y, z][i] {
                            -1 => 0.,
                            1 => 20.,
                            _ => CENTER.0[i],
                        })) + shift;
                        let mut s = s0.clone();
                        for i in 0..2 {
                            s.set_pos(i, s.pos(i) - CENTER + center)
                        }
                        s.wrap_positions();
                        let mut integrator: Box<dyn Integrator> = if half {
                            Box::new(LangevinBAOAB::with_constraints(0., 1., 0, c.clone()))
                        } else {
                            Box::new(VelocityVerlet::with_constraints(c.clone()))
                        };
                        let t = integrator
                            .step_with_stress_options(
                                &mut s,
                                &Zero,
                                &mut NeighborList::new(1.),
                                h,
                                &ModelEvaluation::zero(),
                                ConstraintSolverOptions::TIGHT,
                            )
                            .map_err(err)?
                            .stress
                            .ok_or("no stress")?;
                        work.step(1)?;
                        serde_json::to_writer(&mut *raw,&serde_json::json!({"fixture":id,"image":[x,y,z],"shift":shift,"half_maps":half,"dt":h,"q":s.positions,"v":s.velocities,"stages":t.stages,"interval":t.interval,"endpoint":t.endpoint})).map_err(err)?;
                        writeln!(raw).map_err(err)?;
                        let mut e = Vec3::new(1., 2., 2.) * (1. / 3.);
                        let mut tangent = Vec3::new(2., -1., 0.) * (1. / 5_f64.sqrt());
                        for j in 0..n {
                            let en = e * cs + tangent * a;
                            let tn = tangent * cs - e * a;
                            let (pre, post) = if half { [(1, 2), (4, 5)][j] } else { (0, 1) };
                            for (offset, index, direction) in [(0, pre, e), (1, post, en)] {
                                let expected = dyad(direction, direction)
                                    .map(|v| v * mu * 1.3_f64.powi(2) * (cs - 1.) / tau);
                                let d = maxdiff(&t.stages[index].tensor_impulse.0, &expected)?;
                                errors[j * 2 + offset] = errors[j * 2 + offset].max(d);
                            }
                            e = en;
                            tangent = tn;
                        }
                        for i in 0..2 {
                            let frac = if i == 0 { -m[1] / mass } else { m[0] / mass };
                            let dq = s
                                .box_size
                                .min_image(s.pos(i) - (center + V * h + e * (1.3 * frac)));
                            errors[n * 2] = errors[n * 2].max(maxdiff(&dq.0, &[0.; 3])?);
                            errors[n * 2 + 1] = errors[n * 2 + 1]
                                .max(maxdiff(&s.vel(i).0, &(V + tangent * (1.3 * 0.7 * frac)).0)?);
                        }
                        let trace =
                            t.interval.scalar.ok_or("missing scalar")? * 3. * s.box_size.volume();
                        errors[n * 2 + 2] =
                            errors[n * 2 + 2].max(maxdiff(&[trace], &[mass * V.norm2()])?);
                        if !half && id == "rotor" && [x, y, z] == [0; 3] && shift == Vec3::zero() {
                            let wrong = t.endpoint.kinetic.trace()
                                + t.interval.constraint_average.trace()
                                - mass * V.norm2();
                            g.push(scaled(
                                "vv_rotor_mixed_time_negative_control",
                                &[wrong],
                                &[-0.621075 * (1. - cs) / (1. + cs)],
                                5e-10,
                            )?);
                        }
                    }
                }
            }
        }
        g.push(gate(
            &format!(
                "{id}_{}_maps_all_images",
                if half { "baoab_two_half" } else { "vv" }
            ),
            &errors,
            &vec![0.; errors.len()],
            &vec![5e-10; errors.len()],
        )?);
    }
    Ok(())
}
pub fn run_deterministic(out: &Path) -> R<Vec<GateResult>> {
    use std::io::Write;
    let mut raw = std::io::BufWriter::new(
        std::fs::File::create_new(out.join("deterministic.jsonl")).map_err(err)?,
    );
    let mut gates = vec![];
    let mut work = Work::begin();
    for f in det_fixtures() {
        let s0 = f.system()?;
        let id = &f.id;
        let c = &f.constraints;
        let time = 8. / 0.7;
        if id.starts_with("rotor") {
            finite_maps(&mut gates, &s0, c, id, &mut work, &mut raw)?;
            for tilted in [false, true] {
                let mut trap = s0.clone();
                for i in 0..2 {
                    let mut r = trap.pos(i) - CENTER;
                    let mut u = trap.vel(i) - V;
                    if tilted {
                        r = common::water_rotation(r);
                        u = common::water_rotation(u)
                    }
                    trap.set_pos(i, CENTER + r);
                    trap.set_vel(i, V + u);
                    trap.set_force(i, r * (-0.1));
                }
                let rl = endpoint_reaction(&trap, &constraint_geometry(&trap, c).map_err(err)?)
                    .map_err(err)?;
                let (force, free_and_constraint) = analytic_endpoint(&trap)?;
                let external = modena_core::stress::SymTensor6(std::array::from_fn(|k| {
                    (0..2)
                        .map(|i| dyad(trap.pos(i) - CENTER, trap.force(i))[k])
                        .sum()
                }));
                let m = ModelEvaluation {
                    energy: (0..2).map(|i| 0.05 * (trap.pos(i) - CENTER).norm2()).sum(),
                    scalar_virial: external.trace(),
                    tensor: modena_core::stress::TensorAvailability::Complete(external),
                };
                let pressure =
                    modena_core::stress::endpoint_pressure(&trap, &m, Some(c)).map_err(err)?;
                let volume = trap.box_size.volume();
                let expected = std::array::from_fn::<_, 6, _>(|k| {
                    free_and_constraint[k] * volume + external.0[k]
                });
                gates.push(scaled(
                    &format!("{id}_forced_trap_{tilted}_pressure"),
                    &tensor(&pressure.total)?.map(|p| p * volume),
                    &expected,
                    2e-9,
                )?);
                for (i, (got, want)) in rl.forces.iter().zip(force).enumerate() {
                    gates.push(scaled(
                        &format!("{id}_forced_trap_{tilted}_atom_{i}"),
                        &got.0,
                        &want.0,
                        2e-9,
                    )?);
                }
            }
        }
        let exact = reference(&s0, id, 800 * 64)?;
        if id.contains("nonprincipal") {
            gates.push(gate(
                &format!("{id}_rk4_refinement"),
                &state_err(&reference(&s0, id, 800 * 32)?, &exact)?,
                &[0.; 2],
                &[1e-10; 2],
            )?);
        }
        let g0 = moment(&s0);
        let g1 = moment(&exact);
        let expected: [f64; 6] =
            std::array::from_fn(|k| (g1[k] - g0[k]) / (time * s0.box_size.volume()));
        let mut errors = vec![];
        let mut tensor_errors = vec![];
        for steps in [100, 200, 400, 800] {
            let h = time / steps as f64;
            let mut s = s0.clone();
            let mut integ = VelocityVerlet::with_constraints(c.clone());
            let mut nl = NeighborList::new(1.);
            let mut pi = [0.; 6];
            let mut pe = [0.; 6];
            let mut last = analytic_endpoint(&s)?.1;
            let mut max_f = 0_f64;
            let mut max_p = 0_f64;
            let mut max_res = 0_f64;
            for step in 1..=steps {
                let t = integ
                    .step_with_stress_options(
                        &mut s,
                        &Zero,
                        &mut nl,
                        h,
                        &ModelEvaluation::zero(),
                        ConstraintSolverOptions::TIGHT,
                    )
                    .map_err(err)?
                    .stress
                    .ok_or("missing stress")?;
                work.step(1)?;
                let interval = tensor(&t.interval.total)?;
                let endpoint = tensor(&t.endpoint.total)?;
                let (force, analytic) = analytic_endpoint(&s)?;
                let reaction = endpoint_reaction(&s, &constraint_geometry(&s, c).map_err(err)?)
                    .map_err(err)?;
                for (f, r) in reaction.forces.iter().zip(&force) {
                    for k in 0..3 {
                        let d = (f.0[k] - r.0[k]).abs() / (1. + r.0[k].abs());
                        finite(&[d])?;
                        max_f = max_f.max(d);
                    }
                }
                finite(&[reaction.relative_linear_residual])?;
                max_res = max_res.max(reaction.relative_linear_residual);
                for k in 0..6 {
                    let d = (endpoint[k] - analytic[k]).abs() * s.box_size.volume()
                        / (1. + (analytic[k] * s.box_size.volume()).abs());
                    finite(&[d])?;
                    max_p = max_p.max(d);
                    pi[k] += interval[k] / steps as f64;
                    pe[k] += (endpoint[k] + last[k]) / (2. * steps as f64);
                }
                last = endpoint;
                finite(&s.positions)?;
                finite(&s.velocities)?;
                finite(&pi)?;
                finite(&pe)?;
                serde_json::to_writer(&mut raw,&serde_json::json!({"fixture":id,"h":h,"step":step,"time":step as f64*h,"q":s.positions,"v":s.velocities,"endpoint":endpoint,"interval":interval,"analytic_endpoint":analytic,"reaction_forces":reaction.forces,"analytic_forces":force})).map_err(err)?;
                writeln!(raw).map_err(err)?;
            }
            let st = state_err(&s, &exact)?;
            let ie = maxdiff(&pi, &expected)?;
            let ee = maxdiff(&pe, &expected)?;
            errors.push([st[0], st[1], ie, ee]);
            tensor_errors.push([
                std::array::from_fn::<_, 6, _>(|k| (pi[k] - expected[k]).abs()),
                std::array::from_fn::<_, 6, _>(|k| (pe[k] - expected[k]).abs()),
            ]);
            gates.push(gate(
                &format!("{id}_{steps}_euler_or_rotor_reaction"),
                &[max_f, max_p, max_res],
                &[0.; 3],
                &[2e-9, 2e-9, 1e-11],
            )?);
            let numerical = moment(&s);
            let discrete: [f64; 6] =
                std::array::from_fn(|k| (numerical[k] - g0[k]) / (time * s.box_size.volume()));
            gates.push(scaled(
                &format!("{id}_{steps}_finite_window_ledger"),
                &pi,
                &discrete,
                5e-10 / s.box_size.volume(),
            )?);
            serde_json::to_writer(&mut raw,&serde_json::json!({"fixture":id,"steps":steps,"expected_q":exact.positions,"expected_v":exact.velocities,"expected_window_pressure":expected,"interval_mean":pi,"endpoint_trapezoid_mean":pe,"errors":errors.last()})).map_err(err)?;
            writeln!(raw).map_err(err)?;
        }
        for (kind, label) in ["interval", "endpoint"].iter().enumerate() {
            for component in 0..6 {
                let e = tensor_errors
                    .iter()
                    .map(|p| p[kind][component])
                    .collect::<Vec<_>>();
                finite(&e)?;
                // Resolve individual component order only above a roundoff floor in pressure units.
                if e[3]
                    > 256. * f64::EPSILON * (1. + expected[component].abs()) / s0.box_size.volume()
                {
                    gates.push(gate(
                        &format!("{id}_{label}_{component}_ratios"),
                        &[e[1] / e[2], e[2] / e[3]],
                        &[4.; 2],
                        &[1.; 2],
                    )?);
                }
            }
        }
        for (k, label) in ["positions", "velocities", "interval", "endpoint"]
            .iter()
            .enumerate()
        {
            let ratios = [errors[1][k] / errors[2][k], errors[2][k] / errors[3][k]];
            gates.push(gate(
                &format!("{id}_{label}_last_two_ratios"),
                &ratios,
                &[4.; 2],
                &[1.; 2],
            )?);
        }
    }
    raw.flush().map_err(err)?;
    work.persist(true)?;
    Ok(gates)
}
