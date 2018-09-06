use na::{Matrix2, Vector2};
use num::complex::Complex;

use super::common::{
    rotate_matrix, well_behaved_complexes, well_behaved_doubles, Angle, ComplexMatrix,
    ElementParams, JonesError, JonesMatrix, JonesVector, MissingParameter, Result,
};

#[derive(Debug)]
pub struct Polarizer {
    mat: ComplexMatrix,
}

impl Polarizer {
    pub fn new(angle: Angle) -> Self {
        let rad = match angle {
            Angle::Degrees(ang) => ang.to_radians(),
            Angle::Radians(ang) => ang,
        };
        let cos_2 = Complex::<f64>::new(rad.cos().powi(2), 0.0);
        let sin_2 = Complex::<f64>::new(rad.sin().powi(2), 0.0);
        let sin_cos = Complex::<f64>::new(rad.sin() * rad.cos(), 0.0);
        let mat = Matrix2::new(cos_2, sin_cos, sin_cos, sin_2);
        Polarizer { mat }
    }
}

impl From<ElementParams> for Result<Polarizer> {
    fn from(params: ElementParams) -> Self {
        match params.angle {
            Some(angle) => Ok(Polarizer::new(angle)),
            None => {
                let missing = MissingParameter {
                    typ: "Polarizer".into(),
                    param: "angle".into(),
                };
                Err(JonesError::MissingParameter(missing))
            }
        }
    }
}

impl JonesMatrix for Polarizer {
    fn rotated(&self, angle: Angle) -> ComplexMatrix {
        // Just use the default implementation
        rotate_matrix(&self.matrix(), &angle)
    }

    fn rotate(&mut self, angle: Angle) {
        // Just use the default implementation
        self.mat = rotate_matrix(&self.mat, &angle);
    }

    fn matrix(&self) -> ComplexMatrix {
        self.mat
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_horizontal_polarizer() {
        let pol = Polarizer::new(Angle::Degrees(0 as f64));
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
        let pol = Polarizer::new(Angle::Degrees(90 as f64));
        let expected = ComplexMatrix::new(
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0),
            Complex::new(1.0, 0.0),
        );
        assert_matrix_approx_eq!(pol.matrix(), expected);
    }
}

proptest! {
   #[test]
   fn test_polarizer_attenuation(theta in 0 as f64..90 as f64) {
       let beam = Vector2::new(Complex::<f64>::new(1.0, 0.0), Complex::<f64>::new(0.0, 0.0));
       let pol = Polarizer::new(Angle::Degrees(theta));
       let beam_after_pol = pol.matrix() * beam;
       let expected_intensity = theta.to_radians().cos().powi(2);
       assert_approx_eq!(expected_intensity, beam_after_pol.intensity().unwrap());
   }

   #[test]
   fn test_crossed_polarizers(x in well_behaved_complexes(),
                              y in well_behaved_complexes(),
                              theta in 0 as f64..90 as f64) {
       let beam = Vector2::new(x, y);
       let first_pol = Polarizer::new(Angle::Degrees(theta));
       let second_pol = Polarizer::new(Angle::Degrees(theta + 90.0));
       let beam_after = second_pol.matrix() * first_pol.matrix() * beam;
       assert_approx_eq!(0.0, beam_after.intensity().unwrap());
   }
}
