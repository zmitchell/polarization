//! A linear polarizer.
use na::Matrix2;
use num::complex::Complex;

use super::common::{rotate_matrix, Angle, ComplexMatrix, JonesMatrix};
#[cfg(test)]
use proptest::prelude::*;

/// An ideal linear polarizer.
#[derive(Debug, Copy, Clone)]
pub struct Polarizer {
    mat: ComplexMatrix,
}

impl Polarizer {
    /// Constructs a new linear polarizer oriented at `angle`.
    pub fn new(angle: Angle) -> Self {
        let rad = match angle {
            Angle::Degrees(ang) => ang.to_radians(),
            Angle::Radians(ang) => ang,
        };
        let cos_2 = Complex::new(rad.cos().powi(2), 0_f64);
        let sin_2 = Complex::new(rad.sin().powi(2), 0_f64);
        let sin_cos = Complex::new(rad.sin() * rad.cos(), 0_f64);
        let mat = Matrix2::new(cos_2, sin_cos, sin_cos, sin_2);
        Polarizer { mat }
    }
}

impl JonesMatrix for Polarizer {
    /// Returns the element rotated counter-clockwise by `angle`.
    fn rotated(&self, angle: Angle) -> Self {
        // Just use the default implementation
        Polarizer {
            mat: rotate_matrix(&self.matrix(), &angle),
        }
    }

    /// Rotate the element counter-clockwise by `angle`.
    fn rotate(&mut self, angle: Angle) {
        // Just use the default implementation
        self.mat = rotate_matrix(&self.mat, &angle);
    }

    /// Returns the 2x2 Jones matrix of the element.
    fn matrix(&self) -> ComplexMatrix {
        self.mat
    }
}

#[cfg(test)]
impl Arbitrary for Polarizer {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        any::<Angle>()
            .prop_map(|angle| Polarizer::new(angle))
            .boxed()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use jones::common::{float_angle, Beam, JonesVector};

    #[test]
    fn test_horizontal_polarizer() {
        let pol = Polarizer::new(Angle::Degrees(0_f64));
        let expected = ComplexMatrix::new(
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
        );
        assert_matrix_approx_eq!(pol.matrix(), expected);
    }

    #[test]
    fn test_vertical_polarizer() {
        let pol = Polarizer::new(Angle::Degrees(90_f64));
        let expected = ComplexMatrix::new(
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
        );
        assert_matrix_approx_eq!(pol.matrix(), expected);
    }

    proptest! {
       #[test]
       fn test_polarizer_attenuation(theta in float_angle()) {
           let beam = Beam::new(Complex::new(1_f64, 0_f64), Complex::new(0_f64, 0_f64));
           let pol = Polarizer::new(Angle::Degrees(theta));
           let beam_after = beam.apply_element(pol);
           let expected_intensity = theta.to_radians().cos().powi(2);
           assert_approx_eq!(beam_after.intensity().unwrap(), expected_intensity);
       }

       #[test]
       fn test_crossed_polarizers(beam: Beam, theta in float_angle()) {
           let first_pol = Polarizer::new(Angle::Degrees(theta));
           let second_pol = Polarizer::new(Angle::Degrees(theta + 90.0));
           let beam_after = beam.apply_element(first_pol).apply_element(second_pol);
           assert_approx_eq!(beam_after.intensity().unwrap(), 0.0);
       }

       #[test]
       fn test_horizontal_polarizer_zeroes_y_phase(beam: Beam) {
            let pol = Polarizer::new(Angle::Degrees(0.0));
            let beam_after = beam.apply_element(pol).remove_common_phase();
            let y_phase = beam_after.y().arg();
            prop_assert!(y_phase.abs() < 1e-6);
       }

       #[test]
       fn test_vertical_polarizer_zeroes_x_phase(beam: Beam) {
            let pol = Polarizer::new(Angle::Degrees(90.0));
            let beam_after = beam.apply_element(pol).remove_common_phase();
            let x_phase = beam_after.x().arg();
            prop_assert!(x_phase.abs() < 1e-6);
       }
    }
}
