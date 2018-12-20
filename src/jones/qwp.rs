//! A quarter-wave plate.
//!
//! A quarter-wave plate is an optical retarder commonly used to convert a linearly
//! polarized beam into a circularly polarized beam or vice versa.
use na::Matrix2;
use num::complex::Complex;

use super::common::{rotate_matrix, Angle, ComplexMatrix, JonesMatrix};
#[cfg(test)]
use proptest::prelude::*;

/// An ideal quarter-wave plate.
///
/// See the module-level documentation for more details.
#[derive(Debug, Copy, Clone)]
pub struct QuarterWavePlate {
    mat: ComplexMatrix,
}

impl QuarterWavePlate {
    /// Constructs a new quarter-wave plate with its fast axis at `angle`.
    pub fn new(angle: Angle) -> Self {
        let rad = match angle {
            Angle::Radians(ang) => ang,
            Angle::Degrees(ang) => ang.to_radians(),
        };
        let cos_squared = Complex::new(rad.cos().powi(2), 0_f64);
        let sin_squared = Complex::new(rad.sin().powi(2), 0_f64);
        let sin_cos = Complex::new(rad.sin() * rad.cos(), 0_f64);
        let i = Complex::<f64>::i();
        let mat = Matrix2::new(
            cos_squared + i * sin_squared,
            sin_cos - i * sin_cos,
            sin_cos - i * sin_cos,
            sin_squared + i * cos_squared,
        );
        QuarterWavePlate { mat }
    }
}

impl JonesMatrix for QuarterWavePlate {
    /// Returns the element rotated counter-clockwise by `angle`.
    fn rotated(&self, angle: Angle) -> Self {
        QuarterWavePlate {
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
impl Arbitrary for QuarterWavePlate {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<Angle>()
            .prop_map(|angle| QuarterWavePlate::new(angle))
            .boxed()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use jones::common::{are_rel_eq, pi, Beam, JonesVector};

    #[test]
    fn test_horizontal_qwp_doesnt_modify_horizontal_beam() {
        let beam = Beam::new(Complex::new(1_f64, 0_f64), Complex::new(0_f64, 0_f64));
        let qwp = QuarterWavePlate::new(Angle::Degrees(0.0));
        let beam_after = beam.apply_element(qwp);
        assert_beam_approx_eq!(beam_after, beam);
    }

    #[test]
    fn test_vertical_qwp_doesnt_modify_vertical_beam() {
        let beam = Beam::new(Complex::new(0_f64, 0_f64), Complex::new(0_f64, 0_f64));
        let qwp = QuarterWavePlate::new(Angle::Degrees(90.0));
        let beam_after = beam.apply_element(qwp);
        assert_beam_approx_eq!(beam_after, beam);
    }

    #[test]
    fn test_two_qwps_reflect_polarization() {
        // Horizontal beam, should end up vertical after two QWPs at 45 degrees
        let horiz_beam = Beam::new(Complex::new(1_f64, 0_f64), Complex::new(0_f64, 0_f64));
        let vert_beam = Beam::new(Complex::new(0_f64, 0_f64), Complex::new(1_f64, 0_f64));
        let qwp = QuarterWavePlate::new(Angle::Degrees(45.0));
        let beam_after = horiz_beam.apply_element(qwp).apply_element(qwp);
        assert_beam_approx_eq!(beam_after, vert_beam);
    }

    #[test]
    fn test_qwp_circularly_polarizes_beam() {
        // Horizontal beam, should end up right-hand circularly polarized
        let beam = Beam::new(Complex::new(1_f64, 0_f64), Complex::new(0_f64, 0_f64));
        let qwp = QuarterWavePlate::new(Angle::Degrees(45.0));
        let beam_after = beam.apply_element(qwp);
        // Constant factor out front, e^(i*pi/4) / sqrt(2)
        let prefactor = (Complex::<f64>::i() * pi / 4.0).exp() / (2.0_f64.sqrt());
        // Right-hand circularly polarized beam
        let expected_beam = Beam::new(prefactor, -prefactor * Complex::<f64>::i());
        assert_beam_approx_eq!(beam_after, expected_beam);
    }

    proptest! {
        #[test]
        fn test_qwp_preserves_intensity(qwp: QuarterWavePlate, beam: Beam) {
            let intensity_before = beam.intensity();
            prop_assume!(intensity_before.is_ok());
            let intensity_after = beam.apply_element(qwp).intensity();
            prop_assume!(intensity_after.is_ok());
            prop_assert_approx_eq!(intensity_after.unwrap(), intensity_before.unwrap());
        }
    }
}
