//! Crystal lattice builders.

use crate::system::{SimBox, System};
use crate::vec3::Vec3;

/// Build an FCC lattice with `replicas` unit cells per axis and lattice
/// constant `a`. 4 atoms per cell. Box is orthorhombic `a·replicas`.
pub fn fcc(a: f64, replicas: [u32; 3], element: &str) -> System {
    let basis = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.5, 0.5, 0.0),
        Vec3::new(0.5, 0.0, 0.5),
        Vec3::new(0.0, 0.5, 0.5),
    ];
    let [nx, ny, nz] = replicas;
    let box_size = SimBox::new(Vec3::new(a * nx as f64, a * ny as f64, a * nz as f64));
    let mut system = System::new(box_size);
    let total = (nx * ny * nz * 4) as usize;
    system.resize(total);
    let mut idx = 0;
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                for b in &basis {
                    system.positions[idx] = Vec3::new(
                        a * (ix as f64 + b.0[0]),
                        a * (iy as f64 + b.0[1]),
                        a * (iz as f64 + b.0[2]),
                    );
                    system.elements[idx] = element.to_string();
                    idx += 1;
                }
            }
        }
    }
    system
}
