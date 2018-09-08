use na::{Matrix2, Vector2};
use num::complex::Complex;

use super::common::{
    rotate_matrix, well_behaved_complexes, well_behaved_doubles, Angle, Beam, ComplexMatrix,
    ElementParams, JonesError, JonesMatrix, JonesVector, MissingParameter, Result,
};

#[derive(Debug, Copy, Clone)]
pub struct HalfWavePlate {
    mat: ComplexMatrix,
}

impl HalfWavePlate {
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
    fn rotated(&self, angle: Angle) -> Self {
        HalfWavePlate {
            mat: rotate_matrix(&self.mat, &angle),
        }
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
}

proptest!{
   #[test]
   fn test_hwp_reflects_polarization(theta in 0_f64..90_f64) {
       let beam = Beam::linear(Angle::Degrees(theta));
       let expected_beam = Beam::linear(Angle::Degrees(-theta));
       let hwp = HalfWavePlate::new(Angle::Degrees(0.0));
       let beam_after = beam.apply_element(hwp);
       assert_beam_approx_eq!(beam_after, expected_beam);
   }
}
