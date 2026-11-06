//! Certified, allocation-free periodic geometry; no public graph dependency.
//!
//! Intervals enclose real operations on the supplied binary64 inputs. Rounding
//! envelopes enlarge searches only; the final predicate has no epsilon shell.
use crate::{SimBox, Vec3};

#[derive(Clone, Copy, Debug)]
pub(crate) struct Interval {
    pub(crate) lo: f64,
    pub(crate) hi: f64,
}
#[derive(Clone, Copy, Debug)]
pub(crate) enum GeometryError {
    Invalid,
    Uncertified,
    ShiftRange,
}
#[derive(Clone, Copy)]
pub(crate) struct GeometryOptions {
    pub(crate) cutoff: f64,
    pub(crate) pbc: [bool; 3],
    pub(crate) max_abs_shift: i64,
}
#[derive(Clone, Copy, Debug)]
pub(crate) struct CertifiedCell {
    pub(crate) basis: [Vec3; 3],
    pub(crate) inverse: [Vec3; 3],
    pub(crate) inverse_error_inf: f64,
    pub(crate) inverse_row_norm_upper: [f64; 3],
}

const ZERO: Interval = Interval { lo: 0.0, hi: 0.0 };
const ONE: Interval = Interval { lo: 1.0, hi: 1.0 };
const TINY: f64 = f64::from_bits(1);
const MAX_EXACT_SHIFT: i64 = 1_i64 << 52;

impl Interval {
    fn new(lo: f64, hi: f64) -> Result<Self, GeometryError> {
        if !lo.is_finite() || !hi.is_finite() || lo > hi {
            return Err(GeometryError::Uncertified);
        }
        Ok(Self { lo, hi })
    }
    pub(crate) fn point(x: f64) -> Result<Self, GeometryError> {
        Self::new(x, x)
    }
    pub(crate) fn add(self, rhs: Self) -> Result<Self, GeometryError> {
        Self::new((self.lo + rhs.lo).next_down(), (self.hi + rhs.hi).next_up())
    }
    pub(crate) fn sub(self, rhs: Self) -> Result<Self, GeometryError> {
        Self::new((self.lo - rhs.hi).next_down(), (self.hi - rhs.lo).next_up())
    }
    pub(crate) fn mul(self, rhs: Self) -> Result<Self, GeometryError> {
        let products = [
            self.lo * rhs.lo,
            self.lo * rhs.hi,
            self.hi * rhs.lo,
            self.hi * rhs.hi,
        ];
        Self::endpoint_hull(products)
    }
    pub(crate) fn div(self, rhs: Self) -> Result<Self, GeometryError> {
        if rhs.lo <= 0.0 && rhs.hi >= 0.0 {
            return Err(GeometryError::Uncertified);
        }
        Self::endpoint_hull([
            self.lo / rhs.lo,
            self.lo / rhs.hi,
            self.hi / rhs.lo,
            self.hi / rhs.hi,
        ])
    }
    fn endpoint_hull(values: [f64; 4]) -> Result<Self, GeometryError> {
        if values.iter().any(|v| !v.is_finite()) {
            return Err(GeometryError::Uncertified);
        }
        Self::new(
            values.into_iter().fold(f64::INFINITY, f64::min).next_down(),
            values
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        )
    }
    pub(crate) fn sqrt(self) -> Result<Self, GeometryError> {
        if self.lo < 0.0 {
            return Err(GeometryError::Uncertified);
        }
        // sqrt(0) and underflow must not create a negative-domain lower bound.
        Self::new(
            self.lo.sqrt().next_down().max(0.0),
            self.hi.sqrt().next_up(),
        )
    }
    pub(crate) fn square(self) -> Result<Self, GeometryError> {
        let hi = self.abs_upper();
        let lo = if self.lo <= 0.0 && self.hi >= 0.0 {
            0.0
        } else {
            self.lo.abs().min(self.hi.abs())
        };
        Self::new((lo * lo).next_down().max(0.0), (hi * hi).next_up())
    }
    pub(crate) fn hull(self, rhs: Self) -> Self {
        Self {
            lo: self.lo.min(rhs.lo),
            hi: self.hi.max(rhs.hi),
        }
    }
    pub(crate) fn abs_upper(self) -> f64 {
        self.lo.abs().max(self.hi.abs())
    }
}

fn point(x: f64) -> Result<Interval, GeometryError> {
    Interval::point(x)
}
fn norm(values: [Interval; 3]) -> Result<Interval, GeometryError> {
    let mut sum = ZERO;
    for x in values {
        sum = sum.add(x.square()?)?;
    }
    // A sum of squares is nonnegative even when outward addition at zero
    // introduces a negative subnormal lower endpoint.
    Interval::new(sum.lo.max(0.0), sum.hi)?.sqrt()
}
fn abs_row_sum(row: Vec3) -> Result<Interval, GeometryError> {
    let mut sum = ZERO;
    for x in row.0 {
        sum = sum.add(point(x.abs())?)?;
    }
    Ok(sum)
}
fn validate_position(r: Vec3) -> Result<(), GeometryError> {
    if r.0.iter().any(|x| !x.is_finite()) {
        return Err(GeometryError::Invalid);
    }
    Ok(())
}
fn delta(ri: Vec3, rj: Vec3) -> Result<[Interval; 3], GeometryError> {
    validate_position(ri)?;
    validate_position(rj)?;
    let mut d = [ZERO; 3];
    for (k, x) in d.iter_mut().enumerate() {
        *x = point(rj.0[k])?.sub(point(ri.0[k])?)?;
    }
    Ok(d)
}
fn max_abs(r: Vec3) -> f64 {
    r.0.into_iter().map(f64::abs).fold(0.0, f64::max)
}
fn gamma(m: f64) -> Result<Interval, GeometryError> {
    let numerator = point(m)?.mul(point(f64::EPSILON)?)?;
    // Full interval division uses the LOWER denominator for its upper bound.
    numerator.div(ONE.sub(numerator)?)
}
fn underflow(m: f64) -> Result<Interval, GeometryError> {
    point(m)?.mul(point(TINY)?)
}

pub(crate) fn cutoff_squared(cutoff: f64) -> Result<f64, GeometryError> {
    let square = cutoff * cutoff;
    if !cutoff.is_finite() || cutoff <= 0.0 || !square.is_normal() {
        return Err(GeometryError::Invalid);
    }
    Ok(square)
}
fn fp_radius(cutoff2: f64) -> Result<Interval, GeometryError> {
    // Three products and two additions in the nonnegative squared sum.
    // Gradual underflow contributes at most 5*T; no FTZ arithmetic is used.
    point(cutoff2)?
        .add(underflow(5.0)?)?
        .div(ONE.sub(gamma(5.0)?)?)?
        .sqrt()
}

pub(crate) fn certified_cell(cell: SimBox) -> Result<CertifiedCell, GeometryError> {
    let basis = cell.basis();
    let inverse = cell.inv_basis_rows();
    let volume = cell.volume();
    if basis
        .iter()
        .chain(inverse.iter())
        .flat_map(|v| v.0)
        .any(|x| !x.is_finite())
        || !volume.is_finite()
        || volume <= 0.0
    {
        return Err(GeometryError::Invalid);
    }
    let mut rho: f64 = 0.0;
    let mut a_inf: f64 = 0.0;
    for (i, row) in inverse.iter().enumerate() {
        let mut residual_sum = ZERO;
        for (j, column) in basis.iter().enumerate() {
            let mut product = ZERO;
            for k in 0..3 {
                product = product.add(point(row.0[k])?.mul(point(column.0[k])?)?)?;
            }
            let identity = if i == j { ONE } else { ZERO };
            residual_sum = residual_sum.add(point(identity.sub(product)?.abs_upper())?)?;
        }
        rho = rho.max(residual_sum.hi);
        a_inf = a_inf.max(abs_row_sum(*row)?.hi);
    }
    if rho >= 0.5 {
        return Err(GeometryError::Uncertified);
    }
    // E=I-AH, U=(I-E)^(-1) A. Neumann: ||U-A||inf <=
    // rho/(1-rho) ||A||inf. This also bounds every row's 1-norm error.
    let eta = point(rho)?
        .div(ONE.sub(point(rho)?)?)?
        .mul(point(a_inf)?)?
        .hi;
    let mut row_norm = [0.0; 3];
    for (k, row) in inverse.iter().enumerate() {
        row_norm[k] = norm([point(row.0[0])?, point(row.0[1])?, point(row.0[2])?])?
            .add(point(eta)?)?
            .hi;
    }
    Ok(CertifiedCell {
        basis,
        inverse,
        inverse_error_inf: eta,
        inverse_row_norm_upper: row_norm,
    })
}

pub(crate) fn checked_shift_endpoint(x: f64, max_abs: i64) -> Result<i64, GeometryError> {
    if !x.is_finite()
        || !(0..=MAX_EXACT_SHIFT).contains(&max_abs)
        || x < -(max_abs as f64)
        || x > max_abs as f64
        || x.fract() != 0.0
    {
        return Err(GeometryError::ShiftRange);
    }
    // The guard precedes the cast: neither saturation nor loss of integer
    // precision can occur inside [-2^52,2^52].
    Ok(x as i64)
}

pub(crate) fn global_shift_bound(
    cell: &CertifiedCell,
    ri: Vec3,
    rj: Vec3,
    options: GeometryOptions,
) -> Result<i64, GeometryError> {
    let r2 = cutoff_squared(options.cutoff)?;
    checked_shift_endpoint(0.0, options.max_abs_shift)?;
    validate_position(ri)?;
    validate_position(rj)?;
    if !options.pbc.into_iter().any(|x| x) {
        return Ok(0);
    }
    let d = delta(ri, rj)?;
    let mut u_p: f64 = 0.0;
    for (k, enabled) in options.pbc.into_iter().enumerate() {
        if enabled {
            u_p = u_p.max(
                abs_row_sum(cell.inverse[k])?
                    .add(point(cell.inverse_error_inf)?)?
                    .hi,
            );
        }
    }
    let mut h_p: f64 = 0.0;
    for l in 0..3 {
        let mut sum = ZERO;
        for (k, enabled) in options.pbc.into_iter().enumerate() {
            if enabled {
                sum = sum.add(point(cell.basis[k].0[l].abs())?)?;
            }
        }
        h_p = h_p.max(sum.hi);
    }
    let g10 = gamma(10.0)?;
    let alpha = g10.mul(point(u_p)?)?.mul(point(h_p)?)?;
    if alpha.hi >= 0.5 {
        return Err(GeometryError::Uncertified);
    }
    let p = point(max_abs(ri))?.add(point(max_abs(rj))?)?;
    let d_max = d.into_iter().map(Interval::abs_upper).fold(0.0, f64::max);
    let q0 = fp_radius(r2)?.hi.max(options.cutoff);
    let c = point(u_p)?.mul(
        point(d_max)?
            .add(point(q0)?)?
            .add(g10.mul(p)?)?
            .add(underflow(10.0)?)?,
    )?;
    // Entire i64 domain, not just configured B: |fl(S)-S| <= EPS*|S|
    // (including i64::MIN/MAX). One subtraction, three casts, three products,
    // three additions give ||d-x||inf <= gamma_10*(p+H_P*K)+10*T.
    // Any passing norm is finite, so none of these intermediates was nonfinite.
    // Applying periodic rows of U gives K <= C+alpha*K, hence K<=C/(1-alpha).
    // For ANY real integer lattice shift with ||x||<R (even beyond i64),
    // K<=U_P*(D+R)<=C/(1-alpha). Thus the same bound covers physical images
    // and hypothetical f64-passing i64 shifts with inexact casts, before B
    // enters the proof. Accepted G<=B<=2^52 excludes every unsupported image.
    let bound = c.div(ONE.sub(alpha)?)?.hi.ceil();
    checked_shift_endpoint(bound, options.max_abs_shift)
}

pub(crate) fn image_bounds(
    cell: &CertifiedCell,
    ri: Vec3,
    rj: Vec3,
    options: GeometryOptions,
) -> Result<[(i64, i64); 3], GeometryError> {
    // No local B-dependent envelope is available without this global proof.
    global_shift_bound(cell, ri, rj, options)?;
    let mut bounds = [(0, 0); 3];
    if !options.pbc.into_iter().any(|x| x) {
        return Ok(bounds);
    }
    let d = delta(ri, rj)?;
    let d_max = d.into_iter().map(Interval::abs_upper).fold(0.0, f64::max);
    let inverse_error = point(cell.inverse_error_inf)?.mul(point(d_max)?)?.hi;
    let mut e = [ZERO; 3];
    for (k, error) in e.iter_mut().enumerate() {
        let h_sum = abs_row_sum(Vec3(cell.basis.map(|v| v.0[k])))?;
        let magnitude = point(ri.0[k].abs())?
            .add(point(rj.0[k].abs())?)?
            .add(point(options.max_abs_shift as f64)?.mul(h_sum)?)?;
        // Exactly seven scalar arithmetic operations: subtraction, three
        // products, three additions. The globally certified casts are exact.
        *error = gamma(7.0)?.mul(magnitude)?.add(underflow(7.0)?)?;
    }
    let q = fp_radius(cutoff_squared(options.cutoff)?)?.add(norm(e)?)?;
    for (k, bound) in bounds.iter_mut().enumerate() {
        if !options.pbc[k] {
            continue;
        }
        let mut s = ZERO;
        for (a, x) in cell.inverse[k].0.into_iter().zip(d) {
            s = s.add(point(a)?.mul(x)?)?;
        }
        s = s.add(Interval::new(-inverse_error, inverse_error)?)?;
        let width = q.mul(point(cell.inverse_row_norm_upper[k])?)?;
        let negative_s = ZERO.sub(s)?;
        *bound = (
            checked_shift_endpoint(negative_s.sub(width)?.lo.ceil(), options.max_abs_shift)?,
            checked_shift_endpoint(negative_s.add(width)?.hi.floor(), options.max_abs_shift)?,
        );
    }
    Ok(bounds)
}

pub(crate) fn edge_displacement(
    cell: &CertifiedCell,
    ri: Vec3,
    rj: Vec3,
    shift: [i64; 3],
) -> Result<Vec3, GeometryError> {
    validate_position(ri)?;
    validate_position(rj)?;
    let mut d = [0.0; 3];
    for (k, component) in d.iter_mut().enumerate() {
        let x0 = rj.0[k] - ri.0[k];
        let p0 = cell.basis[0].0[k] * shift[0] as f64;
        let x1 = x0 + p0;
        let p1 = cell.basis[1].0[k] * shift[1] as f64;
        let x2 = x1 + p1;
        let p2 = cell.basis[2].0[k] * shift[2] as f64;
        *component = x2 + p2;
        if [x0, p0, x1, p1, x2, p2, *component]
            .into_iter()
            .any(|x| !x.is_finite())
        {
            return Err(GeometryError::Uncertified);
        }
    }
    Ok(Vec3(d))
}

pub(crate) fn strict_radius(d: Vec3, cutoff2: f64) -> Result<bool, GeometryError> {
    if !cutoff2.is_normal() || cutoff2 <= 0.0 {
        return Err(GeometryError::Invalid);
    }
    validate_position(d)?;
    let n2 = (d.0[0] * d.0[0] + d.0[1] * d.0[1]) + d.0[2] * d.0[2];
    if !n2.is_finite() {
        return Err(GeometryError::Uncertified);
    }
    Ok(n2 < cutoff2)
}

/// Radius for the existing TRICLINIC scalar `SimBox::min_image` expression.
/// The later stencil maps these private errors to its fallback reasons; this
/// helper neither selects a strategy nor changes the scalar physical predicate.
pub(crate) fn scalar_min_image_search_radius(
    cell: &CertifiedCell,
    positions: &[f64],
    radius: f64,
) -> Result<f64, GeometryError> {
    let r2 = cutoff_squared(radius)?;
    let t_error = scalar_fractional_errors(cell, positions)?;
    let coefficient_bound = scalar_coefficient_bounds(cell)?;
    let mut e = [ZERO; 3];
    for (k, error) in e.iter_mut().enumerate() {
        let mut propagated = ZERO;
        let mut magnitude = ZERO;
        for l in 0..3 {
            let h = point(cell.basis[l].0[k].abs())?;
            let coefficient_error = point(f64::EPSILON)?
                .mul(coefficient_bound[l])?
                .add(point(TINY)?)?;
            propagated = propagated.add(h.mul(t_error[l].add(coefficient_error)?)?)?;
            // Coefficients here are the ROUNDED s-n; include their subtraction
            // error in the gamma_5 product/sum magnitude as well.
            magnitude = magnitude.add(h.mul(coefficient_bound[l].add(coefficient_error)?)?)?;
        }
        *error = propagated
            .add(gamma(5.0)?.mul(magnitude)?)?
            .add(underflow(5.0)?)?;
    }
    // If the scalar accepts, its computed vector norm is <=Q_fp. The exact
    // physical image differs by <=||e||2, hence lies inside this search Q.
    Ok(fp_radius(r2)?.add(norm(e)?)?.hi)
}

// Separated to audit the scalar's wrapping domain and inner index bound;
// both helpers remain private and are used by the production radius proof.
fn scalar_fractional_errors(
    cell: &CertifiedCell,
    positions: &[f64],
) -> Result<[Interval; 3], GeometryError> {
    if !positions.len().is_multiple_of(3) || positions.iter().any(|x| !x.is_finite()) {
        return Err(GeometryError::Invalid);
    }
    let mut p = [0.0_f64; 3];
    for r in positions.chunks_exact(3) {
        for k in 0..3 {
            p[k] = p[k].max(r[k].abs());
        }
    }
    let mut d = [ZERO; 3];
    let mut subtraction_error = [ZERO; 3];
    for k in 0..3 {
        d[k] = point(2.0)?.mul(point(p[k])?)?;
        subtraction_error[k] = point(f64::EPSILON)?.mul(d[k])?.add(point(TINY)?)?;
    }
    let d_max = d.into_iter().map(|x| x.hi).fold(0.0, f64::max);
    let mut t_error = [ZERO; 3];
    for (k, error) in t_error.iter_mut().enumerate() {
        let mut propagated = ZERO;
        let mut dot_magnitude = ZERO;
        for l in 0..3 {
            let a = point(cell.inverse[k].0[l].abs())?;
            propagated = propagated.add(a.mul(subtraction_error[l])?)?;
            dot_magnitude = dot_magnitude.add(a.mul(d[l].add(subtraction_error[l])?)?)?;
        }
        // Vec3::dot: three products, two additions. The inverse residual and
        // the coordinate subtraction error are separate propagated terms.
        let dot_error = gamma(5.0)?.mul(dot_magnitude)?.add(underflow(5.0)?)?;
        *error = point(cell.inverse_error_inf)?
            .mul(point(d_max)?)?
            .add(propagated)?
            .add(dot_error)?;
        let computed_t = dot_magnitude.add(dot_error)?;
        if computed_t.add(*error)?.hi >= MAX_EXACT_SHIFT as f64 {
            return Err(GeometryError::ShiftRange);
        }
    }
    Ok(t_error)
}

fn scalar_coefficient_bounds(cell: &CertifiedCell) -> Result<[Interval; 3], GeometryError> {
    // In this domain t-round(t) is exact: Sterbenz for |t|>=.5, identity
    // otherwise. Its coefficient magnitude is <=.5. Interval evaluation
    // encloses each rounded intermediate in the actual initial a*s+b*s+c*s,
    // norm2 and sqrt, including zero/subnormal intermediate domains.
    let half = Interval::new(-0.5, 0.5)?;
    let mut initial = [ZERO; 3];
    for (k, value) in initial.iter_mut().enumerate() {
        *value = point(cell.basis[0].0[k])?
            .mul(half)?
            .add(point(cell.basis[1].0[k])?.mul(half)?)?
            .add(point(cell.basis[2].0[k])?.mul(half)?)?;
    }
    let r_init = norm(initial)?;
    let mut coefficient_bound = [ZERO; 3];
    for (k, coefficient) in coefficient_bound.iter_mut().enumerate() {
        let row = cell.inverse[k].0;
        let a_norm = norm([point(row[0])?, point(row[1])?, point(row[2])?])?;
        let width = a_norm.mul(r_init)?;
        // Existing scalar's 32*EPS margin is reproduced only to bound its
        // enumeration indices. It is not a graph shell or a new tolerance.
        let margin = point(32.0)?
            .mul(point(f64::EPSILON)?)?
            .mul(ONE.add(width)?)?;
        let j = point(0.5)?.add(width)?.add(margin)?.hi.ceil() + 1.0;
        if !j.is_finite() || j > 4096.0 {
            return Err(GeometryError::ShiftRange);
        }
        *coefficient = point(0.5)?.add(point(j)?)?;
    }
    Ok(coefficient_bound)
}

/// Enclose U*r for the exact supplied position, independently of pair wrapping.
/// The inverse residual bounds each row's 1-norm error, hence eta*||r||inf
/// encloses the missing inverse contribution; interval dots cover all five
/// arithmetic roundings and gradual underflow. No floating remainder is used.
pub(crate) fn fractional_position(
    cell: &CertifiedCell,
    r: Vec3,
) -> Result<[Interval; 3], GeometryError> {
    validate_position(r)?;
    let error = point(cell.inverse_error_inf)?.mul(point(max_abs(r))?)?.hi;
    let inflation = point(-error)?.hull(point(error)?);
    let mut result = [ZERO; 3];
    for (k, value) in result.iter_mut().enumerate() {
        let mut dot = ZERO;
        for l in 0..3 {
            dot = dot.add(point(cell.inverse[k].0[l])?.mul(point(r.0[l])?)?)?;
        }
        *value = dot.add(inflation)?;
    }
    Ok(result)
}

/// A closed enclosure must fit entirely inside one half-open wrapped bin.
/// Equal floors certify the integer range reduction, with the integer exactly
/// representable below 2^52. Subtraction/multiplication are outward; any seam
/// or bin-face ambiguity is rejected rather than assigning a rounded midpoint.
pub(crate) fn fractional_bin(s: Interval, m: usize) -> Result<usize, GeometryError> {
    if m == 0 || m > 1_000_000 || s.abs_upper() >= MAX_EXACT_SHIFT as f64 {
        return Err(GeometryError::ShiftRange);
    }
    let integer = s.lo.floor();
    if integer != s.hi.floor() {
        return Err(GeometryError::Uncertified);
    }
    let scaled = s.sub(point(integer)?)?.mul(point(m as f64)?)?;
    let bin = scaled.lo.floor();
    if bin != scaled.hi.floor() || bin < 0.0 || bin >= m as f64 {
        return Err(GeometryError::Uncertified);
    }
    Ok(bin as usize)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn options() -> GeometryOptions {
        GeometryOptions {
            cutoff: 1.5,
            pbc: [true; 3],
            max_abs_shift: 1_000_000,
        }
    }
    fn triangular() -> SimBox {
        SimBox::triclinic([
            Vec3::new(2., 0., 0.),
            Vec3::new(1., 2., 0.),
            Vec3::new(0., 0., 4.),
        ])
    }
    fn contains(i: Interval, x: f64) {
        assert!(i.lo <= x && x <= i.hi, "{i:?} misses {x}");
    }

    #[test]
    fn interval_contains_exact_dyadics() {
        let a = Interval::point(0.5).unwrap();
        let b = Interval::point(0.25).unwrap();
        contains(a.add(b).unwrap(), 0.75);
        contains(a.mul(b).unwrap(), 0.125);
        contains(a.sub(b).unwrap(), 0.25);
        contains(a.div(b).unwrap(), 2.);
        contains(b.sqrt().unwrap(), 0.5);
        let half = Interval::point(f64::from_bits(1)).unwrap().mul(a).unwrap();
        assert!(half.lo <= 0. && half.hi >= f64::from_bits(1));
        let across = a.hull(Interval::point(-0.25).unwrap());
        contains(across.square().unwrap(), 0.);
        contains(across.square().unwrap(), 0.25);
        assert_eq!(across.abs_upper(), 0.5);
        assert!(a.div(across).is_err());
        assert!(Interval::point(-1.).unwrap().sqrt().is_err());
        assert!(Interval::point(f64::INFINITY).is_err());
        assert!(Interval::point(f64::MAX)
            .unwrap()
            .mul(Interval::point(2.).unwrap())
            .is_err());
        contains(
            Interval::point(0.)
                .unwrap()
                .square()
                .unwrap()
                .sqrt()
                .unwrap(),
            0.,
        );
    }
    #[test]
    fn denominator_uses_lower_enclosure() {
        let denominator = Interval::point(1.)
            .unwrap()
            .sub(Interval::point(0.25).unwrap())
            .unwrap();
        let q = Interval::point(1.).unwrap().div(denominator).unwrap();
        assert!(q.hi >= (1. / denominator.lo).next_up());
        assert!(q.hi > 1. / denominator.hi);
    }
    #[test]
    fn identity_inverse_certificate_encloses_zero() {
        let c = certified_cell(SimBox::cubic(1.)).unwrap();
        assert!(c.inverse_error_inf >= 0. && c.inverse_error_inf < 1e-12);
        assert!(c
            .inverse_row_norm_upper
            .into_iter()
            .all(|x| (1. ..1.000000000001).contains(&x)));
    }
    #[test]
    fn skew_inverse_certificate_encloses_triangular_inverse() {
        for perturb in [1.0_f64.next_down(), 1., 1.0_f64.next_up()] {
            let c = certified_cell(SimBox::triclinic([
                Vec3::new(2., 0., 0.),
                Vec3::new(perturb, 2., 0.),
                Vec3::new(0., 0., 4.),
            ]))
            .unwrap();
            let expected = [[0.5, -perturb / 4., 0.], [0., 0.5, 0.], [0., 0., 0.25]];
            for (row, exact) in c.inverse.iter().zip(expected) {
                for (a, u) in row.0.into_iter().zip(exact) {
                    assert!((a - u).abs() <= c.inverse_error_inf);
                }
            }
            assert!(c.inverse_row_norm_upper[0] >= (0.25 + perturb * perturb / 16.).sqrt());
        }
        assert!(certified_cell(SimBox::triclinic([
            Vec3::new(-2., 0., 0.),
            Vec3::new(1., 2., 0.),
            Vec3::new(0., 0., 4.)
        ]))
        .is_ok());
    }
    #[test]
    fn uncertifiable_inverse_fails() {
        assert!(certified_cell(SimBox::cubic(0.)).is_err());
        assert!(certified_cell(SimBox::cubic(f64::NAN)).is_err());
        // Finite, nonsingular basis whose computed residual cannot be enclosed below .5.
        let c = SimBox::triclinic([
            Vec3::new(1., 1., 0.),
            Vec3::new(1., 1.0_f64.next_up(), 0.),
            Vec3::new(0., 0., 1.),
        ]);
        assert!(matches!(certified_cell(c), Err(GeometryError::Uncertified)));
    }
    #[test]
    fn normal_cutoff_required() {
        for r in [0., -1., f64::NAN, f64::INFINITY, f64::MAX, 1e-200, 1e-160] {
            assert!(matches!(cutoff_squared(r), Err(GeometryError::Invalid)));
        }
        assert_eq!(cutoff_squared(2.).unwrap(), 4.);
        // Empty/nonempty scalar inputs must perform the same cutoff validation.
        let c = certified_cell(triangular()).unwrap();
        for p in [&[][..], &[0., 0., 0.][..]] {
            for r in [f64::MAX, 1e-200, 1e-160] {
                assert!(matches!(
                    scalar_min_image_search_radius(&c, p, r),
                    Err(GeometryError::Invalid)
                ));
            }
        }
    }
    #[test]
    fn integer_endpoints_checked_before_cast() {
        let b = 1_i64 << 52;
        for x in [
            f64::NAN,
            f64::INFINITY,
            0.5,
            (b as f64).next_up(),
            i64::MAX as f64,
            i64::MIN as f64,
        ] {
            assert!(matches!(
                checked_shift_endpoint(x, b),
                Err(GeometryError::ShiftRange)
            ));
        }
        assert_eq!(checked_shift_endpoint(b as f64, b).unwrap(), b);
        assert_eq!(checked_shift_endpoint(-(b as f64), b).unwrap(), -b);
        assert!(checked_shift_endpoint(0., -1).is_err());
        assert!(checked_shift_endpoint(0., b + 1).is_err());
    }
    #[test]
    fn global_bound_is_independent_of_configured_b() {
        let c = certified_cell(triangular()).unwrap();
        let ri = Vec3::zero();
        let rj = Vec3::new(4., 1., 0.);
        let o = options();
        let g = global_shift_bound(&c, ri, rj, o).unwrap();
        assert_eq!(
            g,
            global_shift_bound(
                &c,
                ri,
                rj,
                GeometryOptions {
                    max_abs_shift: 2_000_000,
                    ..o
                }
            )
            .unwrap()
        );
        assert!(matches!(
            global_shift_bound(
                &c,
                ri,
                rj,
                GeometryOptions {
                    max_abs_shift: g - 1,
                    ..o
                }
            ),
            Err(GeometryError::ShiftRange)
        ));
        let bounds = image_bounds(&c, ri, rj, o).unwrap();
        let mut passing = 0;
        for x in -16_i64..=16 {
            for y in -16_i64..=16 {
                for z in -16_i64..=16 {
                    // Independent direct Cartesian expression, no production membership helper.
                    let dx = (((rj.0[0] - ri.0[0]) + 2. * x as f64) + y as f64) + 0. * z as f64;
                    let dy =
                        (((rj.0[1] - ri.0[1]) + 0. * x as f64) + 2. * y as f64) + 0. * z as f64;
                    let dz =
                        (((rj.0[2] - ri.0[2]) + 0. * x as f64) + 0. * y as f64) + 4. * z as f64;
                    if (dx * dx + dy * dy) + dz * dz < 2.25 {
                        passing += 1;
                        assert!(x.abs().max(y.abs()).max(z.abs()) <= g);
                        for (s, (lo, hi)) in [x, y, z].into_iter().zip(bounds) {
                            assert!(lo <= s && s <= hi);
                        }
                    }
                }
            }
        }
        assert!(passing > 0);
        for mask in 0..8 {
            let pbc = std::array::from_fn(|k| mask & (1 << k) != 0);
            let masked = GeometryOptions { pbc, ..o };
            let masked_bounds = image_bounds(&c, ri, rj, masked).unwrap();
            for (k, (lo, hi)) in masked_bounds.into_iter().enumerate() {
                if !pbc[k] {
                    assert_eq!((lo, hi), (0, 0));
                }
            }
            let g = global_shift_bound(&c, ri, rj, masked).unwrap();
            for x in -16_i64..=16 {
                for y in -16_i64..=16 {
                    for z in -16_i64..=16 {
                        let shift = [x, y, z];
                        if shift
                            .into_iter()
                            .zip(pbc)
                            .any(|(s, enabled)| !enabled && s != 0)
                        {
                            continue;
                        }
                        let dx = (4. + 2. * x as f64) + y as f64;
                        let dy = 1. + 2. * y as f64;
                        let dz = 4. * z as f64;
                        if (dx * dx + dy * dy) + dz * dz < 2.25 {
                            assert!(x.abs().max(y.abs()).max(z.abs()) <= g);
                            for (s, (lo, hi)) in shift.into_iter().zip(masked_bounds) {
                                assert!(lo <= s && s <= hi);
                            }
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn global_contraction_rejects_uncertifiable_cancellation() {
        // Exact triangular inverse, benign residual, but gamma_10 * ||U|| * ||H|| is enormous.
        let c = certified_cell(SimBox::triclinic([
            Vec3::new(1., 0., 0.),
            Vec3::new(1e8, 1., 0.),
            Vec3::new(0., 0., 1.),
        ]))
        .unwrap();
        assert!(matches!(
            global_shift_bound(&c, Vec3::zero(), Vec3::zero(), options()),
            Err(GeometryError::Uncertified)
        ));
        assert!(matches!(
            image_bounds(&c, Vec3::zero(), Vec3::zero(), options()),
            Err(GeometryError::Uncertified)
        ));
    }
    #[test]
    fn huge_coordinates_reject_before_enumeration() {
        let c = certified_cell(SimBox::cubic(1.)).unwrap();
        let p = Vec3::splat(1e100);
        assert!(matches!(
            image_bounds(&c, p, p, options()),
            Err(GeometryError::ShiftRange) | Err(GeometryError::Uncertified)
        ));
    }
    #[test]
    fn all_nonperiodic_global_shift_bound_is_zero() {
        let c = certified_cell(triangular()).unwrap();
        let o = GeometryOptions {
            pbc: [false; 3],
            max_abs_shift: 0,
            ..options()
        };
        assert_eq!(
            global_shift_bound(&c, Vec3::zero(), Vec3::splat(1e100), o).unwrap(),
            0
        );
        assert_eq!(
            image_bounds(&c, Vec3::zero(), Vec3::splat(1e100), o).unwrap(),
            [(0, 0); 3]
        );
        assert!(global_shift_bound(&c, Vec3::splat(f64::NAN), Vec3::zero(), o).is_err());
        // With no integer axes the shift is identically zero even if a later
        // Cartesian displacement cannot be evaluated finitely.
        assert_eq!(
            global_shift_bound(&c, Vec3::splat(-f64::MAX), Vec3::splat(f64::MAX), o).unwrap(),
            0
        );
        assert_eq!(
            image_bounds(&c, Vec3::splat(-f64::MAX), Vec3::splat(f64::MAX), o).unwrap(),
            [(0, 0); 3]
        );
    }
    #[test]
    fn i64_conversion_error_is_in_global_proof() {
        for n in [
            i64::MIN,
            i64::MIN + 1,
            i64::MAX - 1,
            i64::MAX,
            (1_i64 << 52) - 1,
            (1_i64 << 52) + 1,
            (1_i64 << 53) - 1,
            (1_i64 << 53) + 1,
        ] {
            let error = ((n as f64) as i128 - i128::from(n)).abs() as f64;
            assert!(error <= f64::EPSILON * (n as f64).abs());
        }
        assert_ne!(
            ((1_i64 << 53) + 1) as i128,
            (((1_i64 << 53) + 1) as f64) as i128
        );
        let c = certified_cell(SimBox::cubic(1.)).unwrap();
        // Hand evaluation: p=2^51, D=0, U_P=H_P=1. Ten-rounding
        // gamma contributes just over 5, so ceil(K_global)=ceil(6.5+tiny)=7.
        // An incorrect seven-rounding global envelope instead returns 6.
        let translated = Vec3::splat((1_u64 << 50) as f64);
        assert_eq!(
            global_shift_bound(&c, translated, translated, options()).unwrap(),
            7
        );
        assert!(global_shift_bound(
            &c,
            Vec3::zero(),
            Vec3::new((1_i64 << 53) as f64, 0., 0.),
            GeometryOptions {
                max_abs_shift: 1_i64 << 52,
                ..options()
            }
        )
        .is_err());
    }
    #[test]
    fn cartesian_operation_order_boundary() {
        let c = certified_cell(SimBox::cubic(8.)).unwrap();
        for (x, want) in [
            (2.0_f64.next_down(), true),
            (2., false),
            (2.0_f64.next_up(), false),
        ] {
            let d = edge_displacement(&c, Vec3::zero(), Vec3::new(x, 0., 0.), [0; 3]).unwrap();
            assert_eq!(strict_radius(d, 4.).unwrap(), want);
        }
        for (r, want) in [
            (2.0_f64.next_down(), false),
            (2., false),
            (2.0_f64.next_up(), true),
        ] {
            assert_eq!(
                strict_radius(Vec3::new(2., 0., 0.), cutoff_squared(r).unwrap()).unwrap(),
                want
            );
        }
        let c = certified_cell(SimBox::triclinic([
            Vec3::new(1e16, 0., 0.),
            Vec3::new(-1e16, 2., 0.),
            Vec3::new(1., 0., 4.),
        ]))
        .unwrap();
        assert_eq!(
            edge_displacement(&c, Vec3::zero(), Vec3::new(1., 0., 0.), [1, 1, 1])
                .unwrap()
                .0,
            [1., 2., 4.]
        );
        assert!(strict_radius(Vec3::splat(f64::MAX), 4.).is_err());
        assert!(
            edge_displacement(&c, Vec3::splat(f64::MAX), Vec3::splat(-f64::MAX), [0; 3]).is_err()
        );
    }
    #[test]
    fn scalar_radius_encloses_accepted_skew_pairs_and_rejects_domain() {
        let cell = triangular();
        let c = certified_cell(cell).unwrap();
        let mut p = vec![0., 0., 0.];
        for x in -16..=16 {
            p.extend_from_slice(&[x as f64 / 4., 0.25, 0.5]);
        }
        let q = scalar_min_image_search_radius(&c, &p, 1.5).unwrap();
        assert!((1.5..1.50000001).contains(&q));
        for coordinates in p[3..].chunks_exact(3) {
            let delta = Vec3::new(coordinates[0], coordinates[1], coordinates[2]);
            if cell.min_image(delta).norm2() < 2.25 {
                // Independent fixed cube finds the exact dyadic physical image,
                // without production inverse/wrapping/search-radius helpers.
                let mut physical2 = f64::INFINITY;
                for i in -16..=16 {
                    for j in -16..=16 {
                        for k in -16..=16 {
                            let dx = delta.0[0] + 2. * f64::from(i) + f64::from(j);
                            let dy = delta.0[1] + 2. * f64::from(j);
                            let dz = delta.0[2] + 4. * f64::from(k);
                            physical2 = physical2.min((dx * dx + dy * dy) + dz * dz);
                        }
                    }
                }
                assert!(physical2 <= q * q);
            }
        }
        assert!(scalar_min_image_search_radius(&c, &[1e100; 3], 1.5).is_err());
        assert!(scalar_min_image_search_radius(&c, &[0.; 2], 1.5).is_err());
    }
    #[test]
    fn scalar_t_error_encloses_dyadic_inverse_and_wrap_domain() {
        let c = certified_cell(triangular()).unwrap();
        for scale in [1., 1024., (1_u64 << 40) as f64] {
            let ri = Vec3::new(scale, -scale / 2., scale / 4.);
            let rj = Vec3::new(-scale / 4., scale / 8., -scale / 16.);
            let positions: Vec<_> = ri.0.into_iter().chain(rj.0).collect();
            let errors = scalar_fractional_errors(&c, &positions).unwrap();
            let d = delta(rj, ri).unwrap();
            let exact = [
                d[0].mul(point(0.5).unwrap())
                    .unwrap()
                    .sub(d[1].mul(point(0.25).unwrap()).unwrap())
                    .unwrap(),
                d[1].mul(point(0.5).unwrap()).unwrap(),
                d[2].mul(point(0.25).unwrap()).unwrap(),
            ];
            for k in 0..3 {
                let computed = c.inverse[k].dot(ri - rj);
                assert!(
                    point(computed).unwrap().sub(exact[k]).unwrap().abs_upper() <= errors[k].hi
                );
            }
        }
        let unit = certified_cell(SimBox::cubic(1.)).unwrap();
        assert!(scalar_fractional_errors(&unit, &[(1_u64 << 51) as f64 - 16., 0., 0.]).is_ok());
        assert!(matches!(
            scalar_fractional_errors(&unit, &[(1_u64 << 51) as f64, 0., 0.]),
            Err(GeometryError::ShiftRange)
        ));
        assert!(matches!(
            scalar_fractional_errors(&unit, &[(1_u64 << 52) as f64, 0., 0.]),
            Err(GeometryError::ShiftRange)
        ));
    }

    #[test]
    fn scalar_inner_index_cap_is_certified_before_enumeration() {
        let c = certified_cell(triangular()).unwrap();
        assert!(scalar_coefficient_bounds(&c)
            .unwrap()
            .into_iter()
            .all(|x| x.hi <= 4096.5_f64.next_up()));
        let c = certified_cell(SimBox::triclinic([
            Vec3::new(1., 0., 0.),
            Vec3::new(10_000., 1., 0.),
            Vec3::new(0., 0., 1.),
        ]))
        .unwrap();
        assert!(matches!(
            scalar_coefficient_bounds(&c),
            Err(GeometryError::ShiftRange)
        ));
        assert!(matches!(
            scalar_min_image_search_radius(&c, &[0.; 3], 1.),
            Err(GeometryError::ShiftRange)
        ));
    }

    #[test]
    fn fractional_single_bin_rejects_closed_faces_and_seams() {
        for x in [0.125, -0.875, 2.125] {
            assert_eq!(fractional_bin(point(x).unwrap(), 4).unwrap(), 0);
        }
        for x in [0., 0.25, 1., -1.] {
            assert!(fractional_bin(point(x).unwrap(), 4).is_err());
        }
        assert!(fractional_bin(
            Interval {
                lo: 0.249999999999,
                hi: 0.250000000001
            },
            4
        )
        .is_err());
        assert!(matches!(
            fractional_bin(point((1_u64 << 52) as f64).unwrap(), 4),
            Err(GeometryError::ShiftRange)
        ));
        let c = certified_cell(triangular()).unwrap();
        let r = Vec3::new(0.5, 0.5, 0.5);
        for (s, exact) in fractional_position(&c, r)
            .unwrap()
            .into_iter()
            .zip([0.125, 0.25, 0.125])
        {
            contains(s, exact);
        }
    }

    #[test]
    fn scalar_selected_real_images_fit_q_on_all_fixed_cube_fixtures() {
        #[derive(serde::Deserialize)]
        struct Case {
            cell: [[f64; 3]; 3],
            positions: Vec<[f64; 3]>,
            cutoff: f64,
        }
        let cases: Vec<Case> =
            serde_json::from_str(include_str!("../tests/fixtures/periodic_graph/cases.json"))
                .unwrap();
        let mut accepted = 0;
        for f in cases {
            let cell = SimBox::triclinic(f.cell.map(Vec3));
            let c = certified_cell(cell).unwrap();
            let positions: Vec<_> = f.positions.iter().flatten().copied().collect();
            let q = scalar_min_image_search_radius(&c, &positions, f.cutoff).unwrap();
            for i in 0..f.positions.len() {
                for j in i + 1..f.positions.len() {
                    let ri = Vec3(f.positions[i]);
                    let rj = Vec3(f.positions[j]);
                    let d = ri - rj;
                    let t = c.inverse.map(|row| row.dot(d));
                    let wrap = t.map(f64::round);
                    let s = std::array::from_fn::<_, 3, _>(|k| t[k] - wrap[k]);
                    let [a, b, cc] = c.basis;
                    let mut best = a * s[0] + b * s[1] + cc * s[2];
                    let mut shift = wrap.map(|x| x as i64);
                    let radius = best.norm2().sqrt();
                    // Test-only trace of the unchanged scalar selection records
                    // WHICH image won, not just the independent cube minimum.
                    let bounds = std::array::from_fn::<_, 3, _>(|k| {
                        let width = c.inverse[k].norm() * radius;
                        let margin = 32. * f64::EPSILON * (1. + width);
                        (
                            (s[k] - width - margin).ceil() as i64,
                            (s[k] + width + margin).floor() as i64,
                        )
                    });
                    for x in bounds[0].0..=bounds[0].1 {
                        for y in bounds[1].0..=bounds[1].1 {
                            for z in bounds[2].0..=bounds[2].1 {
                                let candidate = a * (s[0] - x as f64)
                                    + b * (s[1] - y as f64)
                                    + cc * (s[2] - z as f64);
                                if candidate.norm2() < best.norm2() {
                                    best = candidate;
                                    shift = [
                                        wrap[0] as i64 + x,
                                        wrap[1] as i64 + y,
                                        wrap[2] as i64 + z,
                                    ];
                                }
                            }
                        }
                    }
                    assert_eq!(best, cell.min_image(d));
                    if best.norm2() < f.cutoff * f.cutoff {
                        // Every supplied fixed-cube fixture winner lies in
                        // [-16,16]^3. Exact-real Cartesian image via outward
                        // intervals is independent of the scalar inverse path.
                        assert!(shift.into_iter().all(|x| x.abs() <= 16));
                        let mut real = delta(rj, ri).unwrap();
                        for (k, component) in real.iter_mut().enumerate() {
                            for (l, n) in shift.into_iter().enumerate() {
                                *component = component
                                    .sub(
                                        point(f.cell[l][k])
                                            .unwrap()
                                            .mul(point(n as f64).unwrap())
                                            .unwrap(),
                                    )
                                    .unwrap();
                            }
                        }
                        assert!(norm(real).unwrap().hi <= q);
                        accepted += 1;
                    }
                }
            }
        }
        assert!(accepted > 20);
    }
}
