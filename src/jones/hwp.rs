use na::{Matrix2, Vector2};
use num::complex::Complex;

use super::common::{
    beam_lin_pol, rotate_matrix, well_behaved_complexes, well_behaved_doubles, Angle,
    ComplexMatrix, ElementParams, JonesError, JonesMatrix, JonesVector, MissingParameter, Result,
};

#[derive(Debug)]
pub struct HalfWavePlate {
    mat: ComplexMatrix,
}

impl HalfWavePlate {
    pub fn new(angle: Angle) -> Self {
        let rad = match angle {
            Angle::Degrees(deg) => deg.to_radians(),
            Angle::Radians(rad) => rad,
        };
        let sin2 = Complex::<f64>::new((2.0 * rad).sin(), 0.0);
        let cos2 = Complex::<f64>::new((2.0 * rad).cos(), 0.0);
        let mat = Matrix2::new(cos2, sin2, sin2, -cos2);
        HalfWavePlate { mat }
    }
}

impl From<ElementParams> for Result<HalfWavePlate> {
    fn from(params: ElementParams) -> Result<HalfWavePlate> {
        match params.angle {
            Some(angle) => Ok(HalfWavePlate::new(angle)),
            None => {
                let missing = MissingParameter {
                    typ: "HalfWavePlate".into(),
                    param: "angle".into(),
                };
                Err(JonesError::MissingParameter(missing))
            }
        }
    }
}

impl JonesMatrix for HalfWavePlate {
    fn rotated(&self, angle: Angle) -> ComplexMatrix {
        rotate_matrix(&self.mat, &angle)
    }

    fn rotate(&mut self, angle: Angle) {
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
    fn test_hwp_ignores_parallel_beam() {
        let beam = Vector2::new(Complex::<f64>::new(1.0, 0.0), Complex::<f64>::new(0.0, 0.0));
        let hwp = HalfWavePlate::new(Angle::Degrees(0.0));
        let beam_after = hwp.matrix() * beam;
        assert_beam_approx_eq!(beam, beam_after);
    }

    #[test]
    fn test_hwp_ignores_perpendicular_beam() {
        let beam = Vector2::new(Complex::<f64>::new(0.0, 0.0), Complex::<f64>::new(1.0, 0.0));
        let hwp = HalfWavePlate::new(Angle::Degrees(0.0));
        let beam_after = hwp.matrix() * hwp.matrix() * beam;
        assert_beam_approx_eq!(beam, beam_after);
    }
}

proptest!{
   #[test]
   fn test_hwp_reflects_polarization(theta in 0 as f64..90 as f64) {
       let beam = beam_lin_pol(Angle::Degrees(theta));
       let expected_reflection = beam_lin_pol(Angle::Degrees(-theta));
       let hwp = HalfWavePlate::new(Angle::Degrees(0.0));
       let beam_after = hwp.matrix() * beam;
       assert_beam_approx_eq!(expected_reflection, beam_after);
   }
}
