//! Energy-only arbitrary-cell mesh oracle, independent of production SPME.
//!
//! H stores column vectors. The direct DFT is unnormalized, and each even
//! Nyquist axis has independently equiprobable positive/negative lifts. This
//! fixed extension makes G(m)=G(-m) for arbitrary H (negate every lift).
use num_complex::Complex;
use std::f64::consts::PI;

fn count(mesh: [usize; 3]) -> Result<usize, String> {
    if mesh.iter().any(|&n| n < 4) {
        return Err("mesh shape needs >=4 planes".into());
    }
    let n = mesh
        .into_iter()
        .try_fold(1usize, |a, b| a.checked_mul(b))
        .ok_or("mesh work-limit overflow")?;
    if n > 225 {
        return Err("direct DFT work limit: at most 225 points".into());
    }
    Ok(n)
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Independent cofactor inverse: rows of H^-1, not columns.
pub fn inverse(cell: [[f64; 3]; 3]) -> Result<(f64, [[f64; 3]; 3]), String> {
    if cell.iter().flatten().any(|v| !v.is_finite()) {
        return Err("nonfinite cell".into());
    }
    let cof = [
        cross(cell[1], cell[2]),
        cross(cell[2], cell[0]),
        cross(cell[0], cell[1]),
    ];
    let volume = (0..3).map(|a| cell[0][a] * cof[0][a]).sum::<f64>();
    if !volume.is_finite() || volume <= 0. {
        return Err("singular or nonpositive cell".into());
    }
    let inv = cof.map(|v| v.map(|x| x / volume));
    if inv.iter().flatten().any(|v| !v.is_finite()) {
        return Err("nonfinite cell inverse".into());
    }
    Ok((volume, inv))
}
fn coordinates(index: usize, mesh: [usize; 3]) -> [usize; 3] {
    [
        index / (mesh[1] * mesh[2]),
        index / mesh[2] % mesh[1],
        index % mesh[2],
    ]
}
type ComplexModes = Vec<([f64; 2], Complex<f64>)>;

fn spectrum_scaled(
    grid: &[f64],
    cell: [[f64; 3]; 3],
    mesh: [usize; 3],
    kappa: f64,
    lambda: Complex<f64>,
) -> Result<ComplexModes, String> {
    if !lambda.re.is_finite() || !lambda.im.is_finite() || lambda.norm_sqr() == 0. {
        return Err("invalid complex uniform scale".into());
    }
    let n = count(mesh)?;
    if grid.len() != n || grid.iter().any(|v| !v.is_finite()) {
        return Err("grid shape or nonfinite values".into());
    }
    if !kappa.is_finite() || kappa <= 0. {
        return Err("nonfinite/nonpositive kappa".into());
    }
    let (volume, inv) = inverse(cell)?;
    let mut out = Vec::with_capacity(n);
    for m in 0..n {
        let indices = coordinates(m, mesh);
        let mut q = [0.; 2];
        for (r, &charge) in grid.iter().enumerate() {
            let position = coordinates(r, mesh);
            let phase = -2.
                * PI
                * (0..3)
                    .map(|a| indices[a] as f64 * position[a] as f64 / mesh[a] as f64)
                    .sum::<f64>();
            let (sin, cos) = phase.sin_cos();
            q[0] += charge * cos;
            q[1] += charge * sin;
        }
        let mut representatives = vec![[0.; 3]];
        let mut window = 1.;
        for a in 0..3 {
            window *= (4. + 2. * (2. * PI * indices[a] as f64 / mesh[a] as f64).cos()) / 6.;
            let signed = if indices[a] <= mesh[a] / 2 {
                indices[a] as f64
            } else {
                indices[a] as f64 - mesh[a] as f64
            };
            for v in &mut representatives {
                v[a] = signed;
            }
            if mesh[a].is_multiple_of(2) && indices[a] == mesh[a] / 2 {
                let mut negative = representatives.clone();
                for v in &mut negative {
                    v[a] = -signed;
                }
                representatives.extend(negative);
            }
        }
        let mut mean = Complex::new(0., 0.);
        if m != 0 {
            for v in &representatives {
                let k: [f64; 3] = std::array::from_fn(|a| {
                    2. * PI * (0..3).map(|b| inv[b][a] * v[b]).sum::<f64>()
                });
                let k2 = k.iter().map(|x| x * x).sum::<f64>();
                let scaled_k2 = Complex::new(k2, 0.) / (lambda * lambda);
                mean += (-scaled_k2 / (4. * kappa * kappa)).exp() / scaled_k2;
            }
            mean /= representatives.len() as f64;
        }
        out.push((
            q,
            2. * PI / volume / window.powi(2) * mean / (lambda * lambda * lambda),
        ));
    }
    if out.iter().any(|(q, g)| {
        !q[0].is_finite() || !q[1].is_finite() || !g.re.is_finite() || !g.im.is_finite()
    }) {
        return Err("nonfinite DFT/influence".into());
    }
    Ok(out)
}
fn spectrum(
    grid: &[f64],
    cell: [[f64; 3]; 3],
    mesh: [usize; 3],
    kappa: f64,
) -> Result<Vec<([f64; 2], f64)>, String> {
    Ok(
        spectrum_scaled(grid, cell, mesh, kappa, Complex::new(1., 0.))?
            .into_iter()
            .map(|(q, g)| (q, g.re))
            .collect(),
    )
}
pub fn grid_energy(
    grid: &[f64],
    cell: [[f64; 3]; 3],
    mesh: [usize; 3],
    kappa: f64,
) -> Result<f64, String> {
    let u = spectrum(grid, cell, mesh, kappa)?
        .iter()
        .map(|(q, g)| g * (q[0] * q[0] + q[1] * q[1]))
        .sum::<f64>();
    if !u.is_finite() {
        return Err("nonfinite grid energy".into());
    }
    Ok(u)
}
/// Actual complex inverse DFT of GQ, without discarding the imaginary part.
pub fn grid_potential(
    grid: &[f64],
    cell: [[f64; 3]; 3],
    mesh: [usize; 3],
    kappa: f64,
) -> Result<Vec<[f64; 2]>, String> {
    let modes = spectrum(grid, cell, mesh, kappa)?;
    let mut potential = vec![[0.; 2]; modes.len()];
    for (r, p) in potential.iter_mut().enumerate() {
        let position = coordinates(r, mesh);
        for (m, (q, g)) in modes.iter().enumerate() {
            let indices = coordinates(m, mesh);
            let phase = 2.
                * PI
                * (0..3)
                    .map(|a| indices[a] as f64 * position[a] as f64 / mesh[a] as f64)
                    .sum::<f64>();
            let (sin, cos) = phase.sin_cos();
            p[0] += g * (q[0] * cos - q[1] * sin);
            p[1] += g * (q[0] * sin + q[1] * cos);
        }
    }
    if potential.iter().flatten().any(|v| !v.is_finite()) {
        return Err("nonfinite potential".into());
    }
    Ok(potential)
}
fn spread(charges: &[f64], fractional: &[[f64; 3]], mesh: [usize; 3]) -> Result<Vec<f64>, String> {
    let mut grid = vec![0.; count(mesh)?];
    if charges.len() != fractional.len()
        || charges
            .iter()
            .chain(fractional.iter().flatten())
            .any(|v| !v.is_finite())
    {
        return Err("particle shape or nonfinite values".into());
    }
    for (&charge, s) in charges.iter().zip(fractional) {
        let mut weights = [[0.; 4]; 3];
        let mut nodes = [0isize; 3];
        for a in 0..3 {
            let x = s[a].rem_euclid(1.) * mesh[a] as f64;
            let floor = x.floor();
            let t = x - floor;
            nodes[a] = floor as isize;
            weights[a] = [
                (1. - t).powi(3) / 6.,
                (3. * t.powi(3) - 6. * t * t + 4.) / 6.,
                (-3. * t.powi(3) + 3. * t * t + 3. * t + 1.) / 6.,
                t.powi(3) / 6.,
            ];
        }
        for x in 0..4 {
            for y in 0..4 {
                for z in 0..4 {
                    let indices = [x, y, z];
                    let wrapped: [usize; 3] = std::array::from_fn(|a| {
                        (nodes[a] + indices[a] as isize - 1).rem_euclid(mesh[a] as isize) as usize
                    });
                    grid[(wrapped[0] * mesh[1] + wrapped[1]) * mesh[2] + wrapped[2]] +=
                        charge * weights[0][x] * weights[1][y] * weights[2][z];
                }
            }
        }
    }
    Ok(grid)
}
pub fn mesh_energy(
    charges: &[f64],
    fractional: &[[f64; 3]],
    cell: [[f64; 3]; 3],
    mesh: [usize; 3],
    kappa: f64,
) -> Result<f64, String> {
    grid_energy(&spread(charges, fractional, mesh)?, cell, mesh, kappa)
}
/// Energy continuation under H -> lambda H. Only V=lambda³ V0 and
/// k²=k0²/lambda² are analytically continued; Q, modes and window stay fixed.
/// Its imaginary energy gives a subtraction-free independent scalar derivative.
pub fn uniform_mesh_energy(
    charges: &[f64],
    fractional: &[[f64; 3]],
    cell: [[f64; 3]; 3],
    mesh: [usize; 3],
    kappa: f64,
    lambda: Complex<f64>,
) -> Result<Complex<f64>, String> {
    let grid = spread(charges, fractional, mesh)?;
    let u: Complex<f64> = spectrum_scaled(&grid, cell, mesh, kappa, lambda)?
        .iter()
        .map(|(q, g)| g * (q[0] * q[0] + q[1] * q[1]))
        .sum();
    if !u.re.is_finite() || !u.im.is_finite() {
        return Err("nonfinite uniform energy".into());
    }
    Ok(u)
}
pub fn strain_cell(cell: [[f64; 3]; 3], component: usize, eta: f64) -> [[f64; 3]; 3] {
    let (a, b) = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)][component];
    cell.map(|v| {
        let mut out = v;
        if a == b {
            out[a] += eta * v[a];
        } else {
            out[a] += eta * v[b] / 2.;
            out[b] += eta * v[a] / 2.;
        }
        out
    })
}
pub fn energy_strains<F: Fn([[f64; 3]; 3]) -> Result<f64, String>>(
    cell: [[f64; 3]; 3],
    eta: f64,
    f: F,
) -> Result<[f64; 6], String> {
    if !eta.is_finite() || eta <= 0. {
        return Err("nonfinite/nonpositive strain step".into());
    }
    inverse(cell)?;
    let mut out = [0.; 6];
    for (c, x) in out.iter_mut().enumerate() {
        *x = -(f(strain_cell(cell, c, eta))? - f(strain_cell(cell, c, -eta))?) / (2. * eta);
    }
    if out.iter().any(|v| !v.is_finite()) {
        return Err("nonfinite energy derivative".into());
    }
    Ok(out)
}
