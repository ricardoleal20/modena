//! Programmatic configuration: the Rust "modena-language" (fluent builder).
//! See `PLAN_MODENA_PY.md` — one validator serves TOML, Python, and Rust.

mod builder;
pub(crate) mod validation;

pub use builder::SimulationBuilder;

/// Construct indexed nonbonded policy using shortest bond connectivity distance.
/// Explicit exceptions take precedence; conflicting duplicate exceptions fail.
pub fn topology_pair_scales(
    spec: &crate::io::input::BondedTopologySpec,
    n: usize,
) -> crate::Result<Vec<crate::neighbor::PairScale>> {
    use crate::neighbor::{PairScale, PairScaleSet};
    use std::collections::{BTreeMap, HashMap, VecDeque};
    let default = if spec.exclude_bonded {
        [0.0; 3]
    } else {
        [1.0; 3]
    };
    let lj = spec.lj_scales.unwrap_or(default);
    let coulomb = spec.coulomb_scales.unwrap_or(default);
    if lj
        .iter()
        .chain(coulomb.iter())
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(crate::ModenaError::Config(
            "topology scales must be finite and nonnegative".into(),
        ));
    }
    let mut graph = vec![Vec::new(); n];
    for b in &spec.bonds {
        let [i, j] = b.atoms;
        if i >= n || j >= n || i == j {
            return Err(crate::ModenaError::Config(
                "invalid topology bond indices".into(),
            ));
        }
        graph[i].push(j);
        graph[j].push(i);
    }
    let mut pairs = BTreeMap::new();
    for i in 0..n {
        if graph[i].is_empty() {
            continue;
        }
        let mut distances = HashMap::from([(i, 0)]);
        let mut queue = VecDeque::from([i]);
        while let Some(j) = queue.pop_front() {
            if distances[&j] == 3 {
                continue;
            }
            for &k in &graph[j] {
                if !distances.contains_key(&k) {
                    distances.insert(k, distances[&j] + 1);
                    queue.push_back(k);
                    if k > i {
                        let d = distances[&k] - 1;
                        pairs.insert(
                            (i, k),
                            PairScale {
                                i,
                                j: k,
                                lj: lj[d],
                                coulomb: coulomb[d],
                            },
                        );
                    }
                }
            }
        }
    }
    // Preserve legacy angle/dihedral endpoint exclusions when no weighting
    // arrays were supplied. Complete bond graphs additionally derive all paths.
    if spec.exclude_bonded && spec.lj_scales.is_none() && spec.coulomb_scales.is_none() {
        let endpoints = spec
            .angles
            .iter()
            .map(|a| (a.atoms[0], a.atoms[2]))
            .chain(spec.dihedrals.iter().map(|d| (d.atoms[0], d.atoms[3])));
        for (i, j) in endpoints {
            let (i, j) = (i.min(j), i.max(j));
            pairs.entry((i, j)).or_insert(PairScale {
                i,
                j,
                lj: 0.0,
                coulomb: 0.0,
            });
        }
    }
    let explicit: Vec<_> = spec
        .exceptions
        .iter()
        .map(|e| PairScale {
            i: e.atoms[0],
            j: e.atoms[1],
            lj: e.lj,
            coulomb: e.coulomb,
        })
        .collect();
    PairScaleSet::new(&explicit)?.validate(n)?;
    for mut e in explicit {
        if e.i > e.j {
            std::mem::swap(&mut e.i, &mut e.j);
        }
        pairs.insert((e.i, e.j), e);
    }
    let result: Vec<_> = pairs.into_values().collect();
    PairScaleSet::new(&result)?.validate(n)?;
    Ok(result)
}
