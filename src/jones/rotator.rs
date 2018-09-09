use super::common::{
    rotate_matrix, Angle, ComplexMatrix, ElementParams, JonesError, JonesMatrix, MissingParameter,
    Result,
};
use na::Matrix2;
use num::complex::Complex;

#[derive(Debug, Copy, Clone)]
pub struct PolarizationRotator {
    mat: ComplexMatrix,
}

impl PolarizationRotator {
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

impl From<ElementParams> for Result<PolarizationRotator> {
    fn from(params: ElementParams) -> Self {
        match params.angle {
            Some(angle) => Ok(PolarizationRotator::new(angle)),
            None => {
                let missing = MissingParameter {
                    typ: "PolarizationRotator".into(),
                    param: "angle".into(),
                };
                Err(JonesError::MissingParameter(missing))
            }
        }
    }
}

impl JonesMatrix for PolarizationRotator {
    fn rotated(&self, angle: Angle) -> Self {
        // Just use the default implementation
        PolarizationRotator {
            mat: rotate_matrix(&self.matrix(), &angle),
        }
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
mod tests {
    use super::*;
    #[macro_use]
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
