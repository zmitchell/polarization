//! A half-wave plate.
//!
//! A half-wave plate (HWP) is an optical element that delays the component of the beam
//! perpendicular to the "fast" axis of the element by a phase of `pi`, or half of a
//! wavelength. When a linearly polarized beam is incident on a half-wave plate, the
//! polarization of the beam is reflected about the fast axis of the HWP. This type of
//! element may be used as an optical isolator, or as a method of rotating the
//! polarization of a beam without reducing its intensity as would be the case with a
//! polarizer.

use na::Matrix2;
use num::complex::Complex;

use super::common::{rotate_matrix, Angle, ComplexMatrix, JonesMatrix};
#[cfg(test)]
use proptest::prelude::*;

/// An ideal half-wave plate.
///
/// See the module-level documentation for more details.
#[derive(Debug, Copy, Clone)]
pub struct HalfWavePlate {
    mat: ComplexMatrix,
}

impl HalfWavePlate {
    /// Constructs a new half-wave plate with its fast axis at `angle`.
    pub fn new(angle: Angle) -> Self {
        let rad = match angle {
            Angle::Degrees(deg) => deg.to_radians(),
            Angle::Radians(rad) => rad,
        };
        let sin2 = Complex::new((2_f64 * rad).sin(), 0_f64);
        let cos2 = Complex::new((2_f64 * rad).cos(), 0_f64);
        let mat = Matrix2::new(cos2, sin2, sin2, -cos2);
        HalfWavePlate { mat }
    }
}

impl JonesMatrix for HalfWavePlate {
    /// Returns the element rotated counter-clockwise by `angle`.
    fn rotated(&self, angle: Angle) -> Self {
        HalfWavePlate {
            mat: rotate_matrix(&self.mat, &angle),
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
impl Arbitrary for HalfWavePlate {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<Angle>()
            .prop_map(|angle| HalfWavePlate::new(angle))
            .boxed()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use jones::common::{Beam, JonesVector};

    #[test]
    fn test_hwp_ignores_parallel_beam() {
        let beam = Beam::new(Complex::new(1_f64, 0_f64), Complex::new(0_f64, 0_f64));
        let hwp = HalfWavePlate::new(Angle::Degrees(0.0));
        let beam_after = beam.apply_element(hwp);
        assert_beam_approx_eq!(beam_after, beam);
    }

    #[test]
    fn test_hwp_ignores_perpendicular_beam() {
        // This is a vertical beam i.e. (0, 1)
        let beam = Beam::new(Complex::new(0_f64, 0_f64), Complex::new(1_f64, 0_f64));
        let hwp = HalfWavePlate::new(Angle::Degrees(0.0));
        // You have to flip reflect the polarization twice to make Rust happy,
        // even though physically the beam is the same no matter how many times
        // the beam travels through the half-wave plate.
        let beam_after = beam.apply_element(hwp).apply_element(hwp);
        assert_beam_approx_eq!(beam_after, beam);
    }

    proptest!{
        #[test]
        fn test_hwp_reflects_polarization(theta in 0_f64..90_f64) {
            let beam = Beam::linear(Angle::Degrees(theta));
            let expected_beam = Beam::linear(Angle::Degrees(-theta));
            let hwp = HalfWavePlate::new(Angle::Degrees(0.0));
            let beam_after = beam.apply_element(hwp);
            prop_assert_beam_approx_eq!(beam_after, expected_beam);
        }

        #[test]
        fn test_hwp_preserves_intensity(hwp: HalfWavePlate, beam: Beam) {
            let intensity_before = beam.intensity();
            prop_assume!(intensity_before.is_ok());
            let beam_after = beam.apply_element(hwp);
            let intensity_after = beam_after.intensity();
            prop_assume!(intensity_after.is_ok());
            prop_assert_approx_eq!(intensity_after.unwrap(), intensity_before.unwrap());
        }
    }
}
