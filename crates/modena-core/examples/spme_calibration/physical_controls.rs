//! Registered physical recipes and independent auxiliary sums, adapted from
//! spme_calibration integration controls without changing settings/tolerances.
use super::*;
use std::f64::consts::PI;
impl<F: FnMut(PhysicalEvent<'_>) -> Result<(), CalibrationFailure>> Runner<'_, F> {
    fn members(
        &mut self,
        o: &ReferenceOrigin,
        s: &Snapshot,
        p: &ReferenceSettings,
        label: &str,
    ) -> Result<usize, CalibrationFailure> {
        let w = self.aux_admit(o, s, p, label, false)?;
        let mut real = Vec::new();
        let mut reciprocal = Vec::new();
        real.try_reserve_exact(w.real_terms as usize).map_err(io)?;
        reciprocal
            .try_reserve_exact(w.atom_k_terms as usize)
            .map_err(io)?;
        let b = s.cell.map(|l| (p.real_cutoff / l).ceil() as i32 + 1);
        for i in 0..s.positions.len() {
            for j in 0..s.positions.len() {
                for x in -b[0]..=b[0] {
                    for y in -b[1]..=b[1] {
                        for z in -b[2]..=b[2] {
                            self.budget.reference.tick()?;
                            let n = [x, y, z];
                            if i == j && n == [0; 3] {
                                continue;
                            }
                            let r2: f64 = (0..3)
                                .map(|a| {
                                    (s.positions[i][a] - s.positions[j][a]
                                        + f64::from(n[a]) * s.cell[a])
                                        .powi(2)
                                })
                                .sum();
                            if r2 < p.real_cutoff.powi(2) {
                                real.push((i, j, n));
                            }
                        }
                    }
                }
            }
        }
        let b = s
            .cell
            .map(|l| (p.kspace_cutoff * l / (2. * PI)).ceil() as i32);
        for x in -b[0]..=b[0] {
            for y in -b[1]..=b[1] {
                for z in -b[2]..=b[2] {
                    self.budget.reference.tick()?;
                    let m = [x, y, z];
                    let k2: f64 = (0..3)
                        .map(|a| (2. * PI * f64::from(m[a]) / s.cell[a]).powi(2))
                        .sum();
                    if k2 > 0. && k2 < p.kspace_cutoff.powi(2) {
                        reciprocal.push(m);
                    }
                }
            }
        }
        self.store(
            PhysicalRecord {
                label: label.into(),
                origin: o.clone(),
                snapshot: s.clone(),
                settings: p.clone(),
                evaluation: None,
                audit: None,
                membership: Some(Membership { real, reciprocal }),
                usage: self.budget.usage(),
            },
            true,
        )
    }
    fn derivative(
        &mut self,
        o: &ReferenceOrigin,
        s: &Snapshot,
        p: &ReferenceSettings,
        label: &str,
        base: Option<(usize, usize)>,
    ) -> Result<(usize, usize), CalibrationFailure> {
        let members = self.members(o, s, p, &format!("{label}: membership"))?;
        if let Some((expected, _)) = base {
            if self.out.records[members].membership != self.out.records[expected].membership {
                return Err(failed(format!("{label}: cutoff membership changed")));
            }
        }
        let id = self.oracle(o, s, p, label)?;
        let a = self.out.records[id]
            .audit
            .as_ref()
            .ok_or_else(|| invalid("physical audit"))?;
        let m = self.out.records[members]
            .membership
            .as_ref()
            .ok_or_else(|| invalid("physical membership"))?;
        if a.real_included_terms != m.real.len() as u64
            || a.reciprocal_vectors != m.reciprocal.len() as u64
        {
            return Err(failed("independent/oracle membership counts"));
        }
        if let Some((_, previous)) = base {
            let b = self.out.records[previous]
                .audit
                .as_ref()
                .ok_or_else(|| invalid("base audit"))?;
            if a.real_membership_hash != b.real_membership_hash
                || a.reciprocal_membership_hash != b.reciprocal_membership_hash
            {
                return Err(failed("oracle derivative membership hash"));
            }
        }
        Ok((members, id))
    }
    fn half(
        &mut self,
        o: &ReferenceOrigin,
        s: &Snapshot,
        p: &ReferenceSettings,
    ) -> Result<usize, CalibrationFailure> {
        self.aux_admit(o, s, p, "independent half reciprocal", true)?;
        let c = 4. * PI / s.cell.iter().product::<f64>();
        let mut out = Evaluation {
            energy: 0.,
            strain_w: 0.,
            forces: vec![[0.; 3]; s.positions.len()],
        };
        let b = s
            .cell
            .map(|l| (p.kspace_cutoff * l / (2. * PI)).ceil() as i32);
        for x in 0..=b[0] {
            for y in -b[1]..=b[1] {
                for z in -b[2]..=b[2] {
                    self.budget.reference.tick()?;
                    if !(x > 0 || x == 0 && y > 0 || x == 0 && y == 0 && z > 0) {
                        continue;
                    }
                    let k = [
                        2. * PI * f64::from(x) / s.cell[0],
                        2. * PI * f64::from(y) / s.cell[1],
                        2. * PI * f64::from(z) / s.cell[2],
                    ];
                    let k2 = k.iter().map(|v| v * v).sum::<f64>();
                    if k2 >= p.kspace_cutoff.powi(2) {
                        continue;
                    }
                    let pref = c * (-k2 / (4. * p.kappa.powi(2))).exp() / k2;
                    let mut density = s.charges.iter().map(|q| q * q).sum::<f64>();
                    for i in 0..s.positions.len() {
                        for j in 0..i {
                            self.budget.reference.tick()?;
                            let phase = (0..3)
                                .map(|a| k[a] * (s.positions[i][a] - s.positions[j][a]))
                                .sum::<f64>();
                            density += 2. * s.charges[i] * s.charges[j] * phase.cos();
                            let f = 2. * pref * s.charges[i] * s.charges[j] * phase.sin();
                            for (a, ka) in k.iter().enumerate() {
                                out.forces[i][a] += f * ka;
                                out.forces[j][a] -= f * ka;
                            }
                        }
                    }
                    out.energy += pref * density;
                    out.strain_w += pref * density * (1. - k2 / (2. * p.kappa.powi(2)));
                }
            }
        }
        self.store(
            PhysicalRecord {
                label: "independent half reciprocal".into(),
                origin: o.clone(),
                snapshot: s.clone(),
                settings: p.clone(),
                evaluation: Some(out),
                audit: None,
                membership: None,
                usage: self.budget.usage(),
            },
            true,
        )
    }
    fn core(
        &mut self,
        o: &ReferenceOrigin,
        s: &Snapshot,
        p: &ReferenceSettings,
    ) -> Result<usize, CalibrationFailure> {
        use modena_core::{
            EwaldForceModel, ForceModel, NeighborList, PairScale, SimBox, System, Vec3, Virial,
        };
        self.reserve_record()?;
        let mut w = planned_work(s, p)?;
        w.real_terms = 0;
        w.atom_k_terms = mul(w.atom_k_terms, s.positions.len() as u64)?;
        // KVec is four f64; geometric Vec growth is bounded by twice the
        // complete rectangular roster. Other fixed3atom model/list storage is
        // conservatively reserved separately; global allocator covers opaque heap.
        let modes = w.atom_k_terms / (2 * (s.positions.len() as u64).pow(2));
        self.heap(add(mul(modes, 64)?, 65536)?)?;
        self.admit(o, "secondary core Ewald", &w, 0, true)?;
        self.budget.accuracy_calls = add(self.budget.accuracy_calls, 1)?;
        self.budget.candidate_real_pair_upper_bound =
            add(self.budget.candidate_real_pair_upper_bound, 3)?;
        self.budget.reference.work.atom_k_terms =
            add(self.budget.reference.work.atom_k_terms, w.atom_k_terms)?;
        let mut system = System::new(SimBox::from_array(s.cell));
        system.resize(s.positions.len());
        for (i, r) in s.positions.iter().enumerate() {
            system.set_pos(i, Vec3(*r));
        }
        system.charges = s.charges.clone();
        let scales: Vec<_> = s
            .scales
            .iter()
            .map(|c| PairScale {
                i: c.i,
                j: c.j,
                lj: 1.,
                coulomb: c.coulomb,
            })
            .collect();
        let model = EwaldForceModel::new(p.kappa, p.real_cutoff, p.kspace_cutoff)
            .map_err(io)?
            .with_pair_scales(scales.clone())
            .map_err(io)?;
        let mut neighbors =
            NeighborList::with_pair_scales(p.real_cutoff, 0.1, &scales).map_err(io)?;
        let mut virial = Virial::default();
        system.zero_forces();
        // Production core has no callback: bounded3atom call, checkpoints outside.
        self.budget.reference.checkpoint()?;
        let energy = model
            .compute(&mut system, &mut virial, &mut neighbors)
            .map_err(io)?;
        let value = Evaluation {
            energy,
            strain_w: virial.w,
            forces: system
                .forces
                .chunks_exact(3)
                .map(|f| [f[0], f[1], f[2]])
                .collect(),
        };
        self.store(
            PhysicalRecord {
                label: "secondary core Ewald".into(),
                origin: o.clone(),
                snapshot: s.clone(),
                settings: p.clone(),
                evaluation: Some(value),
                audit: None,
                membership: None,
                usage: self.budget.usage(),
            },
            false,
        )
    }
}
pub(super) fn special<F: FnMut(PhysicalEvent<'_>) -> Result<(), CalibrationFailure>>(
    r: &mut Runner<'_, F>,
    o: &ReferenceOrigin,
    s: &Snapshot,
) -> Result<(), CalibrationFailure> {
    let p = settings();
    let base = r.derivative(o, s, &p, "force central", None)?;
    let central = r.value(base.1)?.clone();
    for h in [1e-5 * A, 0.5e-5 * A] {
        for atom in 0..3 {
            for axis in 0..3 {
                let mut plus = s.clone();
                plus.positions[atom][axis] += h;
                let mut minus = s.clone();
                minus.positions[atom][axis] -= h;
                let label = format!("force atom{atom} axis{axis} h{h:.17e}");
                let ep = r.derivative(o, &plus, &p, &format!("{label} plus"), Some(base))?;
                let em = r.derivative(o, &minus, &p, &format!("{label} minus"), Some(base))?;
                let fd = -(r.value(ep.1)?.energy - r.value(em.1)?.energy) / (2. * h);
                let analytic = central.forces[atom][axis];
                r.check(
                    label,
                    analytic,
                    fd,
                    1e-7 / A.powi(2) + 1e-6 * analytic.abs().max(fd.abs()),
                )?;
            }
        }
    }
    let base = r.derivative(o, s, &p, "strain central", None)?;
    let central = r.value(base.1)?.clone();
    for h in [1e-5, 0.5e-5] {
        let strain = |lambda: f64| {
            let mut t = s.clone();
            t.cell = t.cell.map(|x| x * lambda);
            for x in &mut t.positions {
                *x = x.map(|v| v * lambda);
            }
            t
        };
        let ep = r.derivative(
            o,
            &strain(1. + h),
            &p,
            &format!("strain plus h{h:.17e}"),
            Some(base),
        )?;
        let em = r.derivative(
            o,
            &strain(1. - h),
            &p,
            &format!("strain minus h{h:.17e}"),
            Some(base),
        )?;
        let fd = -(r.value(ep.1)?.energy - r.value(em.1)?.energy) / (2. * h);
        r.check(
            format!("strain h{h:.17e}"),
            central.strain_w,
            fd,
            1e-7 / A + 1e-6 * central.strain_w.abs().max(fd.abs()),
        )?;
    }
    let full = r.oracle(o, s, &p, "full positive negative k")?;
    let base = r.oracle(
        o,
        s,
        &ReferenceSettings {
            kspace_cutoff: 0.1,
            ..p.clone()
        },
        "zero reciprocal modes",
    )?;
    let half = r.half(o, s, &p)?;
    let mut expected = r.value(half)?.clone();
    let base = r.value(base)?;
    expected.energy += base.energy;
    expected.strain_w += base.strain_w;
    for (f, b) in expected.forces.iter_mut().zip(&base.forces) {
        for a in 0..3 {
            f[a] += b[a];
        }
    }
    r.compare("full half k", &r.value(full)?.clone(), &expected, 1e-10)?;
    let bg = r.oracle(o, s, &p, "charged background")?;
    let audit = r.out.records[bg]
        .audit
        .as_ref()
        .ok_or_else(|| invalid("background audit"))?
        .clone();
    let q = s.charges.iter().sum::<f64>();
    let expected = -PI * q * q / (2. * s.cell.iter().product::<f64>() * p.kappa * p.kappa);
    r.check(
        "background coefficient",
        audit.background_energy,
        expected,
        1e-13 / A,
    )?;
    r.check(
        "background W=3U",
        audit.background_w,
        3. * expected,
        1e-13 / A,
    )?;
    if s.scales.is_empty() {
        let p = ReferenceSettings {
            kappa: 0.5,
            real_cutoff: 3.7,
            kspace_cutoff: 3.71,
        };
        let base = r.oracle(o, s, &p, "bare fraction base")?;
        let base = r.value(base)?.clone();
        for j in [1, 2] {
            for c in [0., 0.4, 0.7, 1.] {
                let mut scaled = s.clone();
                scaled.scales = vec![Scale {
                    i: 0,
                    j,
                    coulomb: c,
                }];
                let got = r.oracle(o, &scaled, &p, &format!("bare fraction j{j} c{c}"))?;
                let d: [f64; 3] = std::array::from_fn(|a| {
                    let v = s.positions[0][a] - s.positions[j][a];
                    v - s.cell[a] * (v / s.cell[a]).round()
                });
                let distance = d.iter().map(|x| x * x).sum::<f64>().sqrt();
                if (distance < p.real_cutoff) != (j == 1) {
                    return Err(failed("bare pair cutoff geometry"));
                }
                let removed = (1. - c) * s.charges[0] * s.charges[j] / distance;
                let mut expected = base.clone();
                expected.energy -= removed;
                expected.strain_w -= removed;
                for (a, da) in d.iter().enumerate() {
                    let f = removed * da / (distance * distance);
                    expected.forces[0][a] -= f;
                    expected.forces[j][a] += f;
                }
                r.compare("bare fraction", &r.value(got)?.clone(), &expected, 1e-11)?;
            }
        }
    }
    let base = r.oracle(o, s, &p, "rigid translation base")?;
    let base = r.value(base)?.clone();
    for shift in r.p.translations.clone() {
        let mut moved = s.clone();
        for x in &mut moved.positions {
            for a in 0..3 {
                x[a] += shift[a] * s.cell[a] / 16.;
            }
        }
        let got = r.oracle(o, &moved, &p, &format!("rigid shift {shift:?}"))?;
        r.compare("rigid translation", &base, &r.value(got)?.clone(), 1e-9)?;
    }
    let base = r.oracle(o, s, &p, "integer wrap base")?;
    let base = r.value(base)?.clone();
    for atom in 0..3 {
        for a in 0..3 {
            for sign in [-1., 1.] {
                let mut moved = s.clone();
                moved.positions[atom][a] += sign * s.cell[a];
                let got = r.oracle(
                    o,
                    &moved,
                    &p,
                    &format!("integer wrap atom{atom} axis{a} sign{sign}"),
                )?;
                r.compare("integer wrap", &base, &r.value(got)?.clone(), 1e-9)?;
            }
        }
    }
    let mut levels = Vec::new();
    for (real_cutoff, kspace_cutoff) in [(8.8, 22.), (11., 24.), (13.2, 26.)] {
        levels.push(r.oracle(
            o,
            s,
            &ReferenceSettings {
                kappa: 2.,
                real_cutoff,
                kspace_cutoff,
            },
            "secondary ladder",
        )?);
    }
    let final_value = r.value(levels[2])?.clone();
    r.strict(
        "secondary final refinement",
        &r.value(levels[1])?.clone(),
        &final_value,
        &final_value,
    )?;
    r.strict(
        "secondary earlier refinement",
        &r.value(levels[0])?.clone(),
        &r.value(levels[1])?.clone(),
        &final_value,
    )?;
    for kappa in [1.8, 2.2] {
        let got = r.oracle(
            o,
            s,
            &ReferenceSettings {
                kappa,
                real_cutoff: 13.2,
                kspace_cutoff: 26.,
            },
            "secondary cross kappa",
        )?;
        r.strict(
            "secondary cross kappa",
            &r.value(got)?.clone(),
            &final_value,
            &final_value,
        )?;
    }
    let short = ReferenceSettings {
        kappa: 2.,
        real_cutoff: 4.41,
        kspace_cutoff: 26.,
    };
    let got = r.oracle(o, s, &short, "secondary omitted image tail")?;
    r.strict(
        "secondary omitted image tail",
        &r.value(got)?.clone(),
        &final_value,
        &final_value,
    )?;
    let core = r.core(o, s, &short)?;
    r.strict(
        "secondary core Ewald",
        &r.value(core)?.clone(),
        &final_value,
        &final_value,
    )
}
pub(super) fn self_images<F: FnMut(PhysicalEvent<'_>) -> Result<(), CalibrationFailure>>(
    r: &mut Runner<'_, F>,
    o: &ReferenceOrigin,
    original: &Snapshot,
) -> Result<(), CalibrationFailure> {
    let mut s = original.clone();
    s.cell = [A; 3];
    s.positions = vec![[0.4, 0.5, 0.6]];
    s.charges = vec![1.];
    s.scales.clear();
    let p = ReferenceSettings {
        kappa: 0.5,
        real_cutoff: 2.42,
        kspace_cutoff: 0.1,
    };
    let got = r.oracle(o, &s, &p, "self images")?;
    let e = r.value(got)?.clone();
    let u_image = 3. * libm::erfc(p.kappa * A) / A;
    let w_image = 3.
        * (libm::erfc(p.kappa * A) / A + 2. * p.kappa / PI.sqrt() * (-(p.kappa * A).powi(2)).exp());
    let bg = -PI / (2. * A.powi(3) * p.kappa.powi(2));
    r.check(
        "six periodic self energies",
        e.energy,
        u_image - p.kappa / PI.sqrt() + bg,
        1e-13 / A,
    )?;
    r.check(
        "six periodic self W",
        e.strain_w,
        w_image + 3. * bg,
        1e-13 / A,
    )?;
    for f in e.forces[0] {
        r.check("self image force cancellation", f, 0., 1e-13 / A.powi(2))?;
    }
    Ok(())
}
