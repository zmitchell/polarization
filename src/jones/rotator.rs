//! A polarization rotator.
//!
//! This optical element may represent a [Faraday rotator](https://en.wikipedia.org/wiki/Faraday_rotator)
//! for example.

use super::common::{rotate_matrix, Angle, ComplexMatrix, JonesMatrix};
use na::Matrix2;
use num::complex::Complex;

/// An optical element that rotates the polarization of a beam.
///
/// The intensity of the beam is preserved by the rotation. See the module-level
/// documentation for more details.
#[derive(Debug, Copy, Clone)]
pub struct PolarizationRotator {
    mat: ComplexMatrix,
}

impl PolarizationRotator {
    /// Constructs a new `PolarizationRotator` that will rotate a beam counter-clockwise
    /// by `angle`.
    pub fn new(angle: Angle) -> Self {
        let rad = match angle {
            Angle::Degrees(ang) => ang.to_radians(),
            Angle::Radians(ang) => ang,
        };
        let sin = Complex::new(rad.sin(), 0_f64);
        let cos = Complex::new(rad.cos(), 0_f64);
        let mat = Matrix2::new(cos, -sin, sin, cos);
        PolarizationRotator { mat }
    }
}

impl JonesMatrix for PolarizationRotator {
    /// Returns the element rotated counter-clockwise by `angle`.
    fn rotated(&self, angle: Angle) -> Self {
        PolarizationRotator {
            mat: rotate_matrix(&self.matrix(), &angle),
        }
    }

    /// Rotate the element counter-clockwise by `angle`.
    fn rotate(&mut self, angle: Angle) {
        self.mat = rotate_matrix(&self.mat, &angle);
    }

    /// Returns the 2x2 Jones matrix of the element.
    fn matrix(&self) -> ComplexMatrix {
        self.mat
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use jones::common::{Beam, JonesVector};

    proptest! {
        #[test]
        fn test_polarization_rotator_rotates(theta in 0 as f64..360 as f64) {
           let beam = Beam::linear(Angle::Degrees(0.0));
           let expected_beam = Beam::linear(Angle::Degrees(theta));
           let rotator = PolarizationRotator::new(Angle::Degrees(theta));
           let beam_after = beam.apply_element(rotator);
           assert_beam_approx_eq!(beam_after, expected_beam);
        }
    }
}
