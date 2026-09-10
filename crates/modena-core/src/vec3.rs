//! Minimal 3-vector type for MD.
//!
//! Stored as array-of-structs (`[f64; 3]`) for clarity in this PoC. The
//! roadmap lists a future SoA + SIMD (`wide`/portable f64x4) layout pass for
//! the hot force kernels; the public `Vec3` API is intentionally small so that
//! migration is mechanical.

use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Vec3(pub [f64; 3]);

impl Vec3 {
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self([x, y, z])
    }

    #[inline]
    pub const fn zero() -> Self {
        Self([0.0, 0.0, 0.0])
    }

    #[inline]
    pub fn splat(v: f64) -> Self {
        Self([v, v, v])
    }

    #[inline]
    pub fn dot(self, other: Self) -> f64 {
        self.0[0] * other.0[0] + self.0[1] * other.0[1] + self.0[2] * other.0[2]
    }

    #[inline]
    pub fn norm2(self) -> f64 {
        self.dot(self)
    }

    #[inline]
    pub fn norm(self) -> f64 {
        self.norm2().sqrt()
    }

    /// Apply minimum-image wrapping for an orthorhombic box along each axis:
    /// `dx -= box * round(dx / box)`.
    #[inline]
    pub fn min_image(self, box_size: Self) -> Self {
        let mut r = self.0;
        for (r_k, b_k) in r.iter_mut().zip(box_size.0.iter()) {
            *r_k -= b_k * (*r_k / b_k).round();
        }
        Self(r)
    }

    #[inline]
    pub fn into_inner(self) -> [f64; 3] {
        self.0
    }
}

impl Add for Vec3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self([
            self.0[0] + rhs.0[0],
            self.0[1] + rhs.0[1],
            self.0[2] + rhs.0[2],
        ])
    }
}

impl Sub for Vec3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self([
            self.0[0] - rhs.0[0],
            self.0[1] - rhs.0[1],
            self.0[2] - rhs.0[2],
        ])
    }
}

impl Neg for Vec3 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self([-self.0[0], -self.0[1], -self.0[2]])
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn mul(self, s: f64) -> Self {
        Self([self.0[0] * s, self.0[1] * s, self.0[2] * s])
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.0[0] += rhs.0[0];
        self.0[1] += rhs.0[1];
        self.0[2] += rhs.0[2];
    }
}

impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        self.0[0] -= rhs.0[0];
        self.0[1] -= rhs.0[1];
        self.0[2] -= rhs.0[2];
    }
}

impl MulAssign<f64> for Vec3 {
    #[inline]
    fn mul_assign(&mut self, s: f64) {
        self.0[0] *= s;
        self.0[1] *= s;
        self.0[2] *= s;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn min_image_wraps() {
        let b = Vec3::new(10.0, 10.0, 10.0);
        assert![(Vec3::new(6.0, 0.0, 0.0).min_image(b) - Vec3::new(-4.0, 0.0, 0.0)).norm() < 1e-12];
        assert![(Vec3::new(-6.0, 0.0, 0.0).min_image(b) - Vec3::new(4.0, 0.0, 0.0)).norm() < 1e-12];
        assert![(Vec3::new(3.0, 3.0, 3.0).min_image(b) - Vec3::new(3.0, 3.0, 3.0)).norm() < 1e-12];
    }
}
