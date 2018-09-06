use na::Matrix2;
use num::complex::Complex;
#[macro_use]
use super::common::{
    ComplexMatrix, ElementParams, JonesError, JonesMatrix, MissingParameter, Result, rotate_matrix,
    beam_lin_pol, Angle,
};

#[derive(Debug)]
pub struct PolarizationRotator {
    mat: ComplexMatrix,
}

impl PolarizationRotator {
    pub fn new(angle: Angle) -> Self {
        let rad = match angle {
            Angle::Degrees(ang) => ang.to_radians(),
            Angle::Radians(ang) => ang,
        };
        let sin = Complex::<f64>::new(rad.sin(), 0.0);
        let cos = Complex::<f64>::new(rad.cos(), 0.0);
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

proptest! {
   #[test]
   fn test_polarization_rotator_rotates(theta in 0 as f64..360 as f64) {
       let beam = beam_lin_pol(Angle::Degrees(0.0));
       let expected_beam = beam_lin_pol(Angle::Degrees(theta));
       let rotator = PolarizationRotator::new(Angle::Degrees(theta));
       let beam_after = rotator.matrix() * beam;
       assert_beam_approx_eq!(beam_after, expected_beam);
   }
}
