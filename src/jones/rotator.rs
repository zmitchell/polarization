//! A polarization rotator.
//!
//! This optical element may represent a [Faraday rotator](https://en.wikipedia.org/wiki/Faraday_rotator)
//! for example.

use super::common::{rotate_matrix, Angle, ComplexMatrix, JonesMatrix};
use na::Matrix2;
use num::complex::Complex;
#[cfg(test)]
use proptest::prelude::*;

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
impl Arbitrary for PolarizationRotator {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<Angle>()
            .prop_map(|angle| PolarizationRotator::new(angle))
            .boxed()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use jones::common::{are_rel_eq, Beam, JonesVector};

    proptest! {
        #[test]
        fn test_polarization_rotator_rotates(theta in 0_f64..360_f64) {
           let beam = Beam::linear(Angle::Degrees(0.0));
           let expected_beam = Beam::linear(Angle::Degrees(theta));
           let rotator = PolarizationRotator::new(Angle::Degrees(theta));
           let beam_after = beam.apply_element(rotator);
           prop_assert_beam_approx_eq!(beam_after, expected_beam);
        }

        #[test]
        fn test_polarization_rotator_preserves_beam_with_0_degree_rotation(beam: Beam) {
            let rotator = PolarizationRotator::new(Angle::Degrees(0.0));
            let beam_after = beam.apply_element(rotator);
            prop_assert_beam_approx_eq!(beam_after, beam);
        }

        #[test]
        fn test_polarization_rotator_preserves_beam_with_360_degree_rotation(beam: Beam) {
            let rotator = PolarizationRotator::new(Angle::Degrees(360.0));
            let beam_after = beam.apply_element(rotator);
            prop_assert_beam_approx_eq!(beam_after, beam);
        }
    }
}
