//! A linear polarizer.
use na::Matrix2;
use num::complex::Complex;

use super::common::{
    rotate_matrix, Angle, ComplexMatrix, JonesMatrix,
};

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
mod test {
    use super::*;
    use jones::common::{well_behaved_complexes, Beam, JonesVector};

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
       fn test_polarizer_attenuation(theta in 0_f64..90_f64) {
           let beam = Beam::new(Complex::new(1_f64, 0_f64), Complex::new(0_f64, 0_f64));
           let pol = Polarizer::new(Angle::Degrees(theta));
           let beam_after = beam.apply_element(pol);
           let expected_intensity = theta.to_radians().cos().powi(2);
           assert_approx_eq!(beam_after.intensity().unwrap(), expected_intensity);
       }

       #[test]
       fn test_crossed_polarizers(x in well_behaved_complexes(),
                                  y in well_behaved_complexes(),
                                  theta in 0_f64..90_f64) {
           let beam = Beam::new(x, y);
           let first_pol = Polarizer::new(Angle::Degrees(theta));
           let second_pol = Polarizer::new(Angle::Degrees(theta + 90.0));
           let beam_after = beam.apply_element(first_pol).apply_element(second_pol);
           assert_approx_eq!(beam_after.intensity().unwrap(), 0.0);
       }
    }
}
