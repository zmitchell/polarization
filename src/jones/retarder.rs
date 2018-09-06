use na::Matrix2;
use num::complex::Complex;

use super::common::{
    beam_lin_pol, rotate_matrix, Angle, ComplexMatrix, ElementParams, JonesError, JonesMatrix,
    MissingParameter, Result,
};

#[derive(Debug)]
pub struct Retarder {
    mat: ComplexMatrix,
}

impl Retarder {
    pub fn new(angle: Angle, phase: Angle) -> Self {
        let angle_rad = match angle {
            Angle::Degrees(deg) => deg.to_radians(),
            Angle::Radians(rad) => rad,
        };
        let phase_rad = match phase {
            Angle::Degrees(deg) => deg.to_radians(),
            Angle::Radians(rad) => rad,
        };
        let sin = Complex::<f64>::new(angle_rad.sin(), 0.0);
        let cos = Complex::<f64>::new(angle_rad.cos(), 0.0);
        let sin_2 = Complex::<f64>::new(angle_rad.sin().powi(2), 0.0);
        let cos_2 = Complex::<f64>::new(angle_rad.cos().powi(2), 0.0);
        let xi = (Complex::<f64>::i() * phase_rad).exp();
        let a = cos_2 + xi * sin_2;
        let b = sin * cos - xi * sin * cos;
        let c = sin * cos - xi * sin * cos;
        let d = sin_2 + xi * cos_2;
        let mat = Matrix2::new(a, b, c, d);
        Retarder { mat }
    }
}

impl From<ElementParams> for Result<Retarder> {
    fn from(params: ElementParams) -> Result<Retarder> {
        let orientation = match params.angle {
            Some(angle) => match angle {
                Angle::Degrees(deg) => Ok(deg.to_radians()),
                Angle::Radians(rad) => Ok(rad),
            },
            None => {
                let missing = MissingParameter {
                    typ: "Retarder".into(),
                    param: "angle".into(),
                };
                Err(JonesError::MissingParameter(missing))
            }
        }?;
        let phase = match params.phase {
            Some(angle) => match angle {
                Angle::Degrees(deg) => Ok(deg.to_radians()),
                Angle::Radians(rad) => Ok(rad),
            },
            None => {
                let missing = MissingParameter {
                    typ: "Retarder".into(),
                    param: "phase".into(),
                };
                Err(JonesError::MissingParameter(missing))
            }
        }?;
        Ok(Retarder::new(
            Angle::Radians(orientation),
            Angle::Radians(phase),
        ))
    }
}

impl JonesMatrix for Retarder {
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

proptest! {
    #[test]
    fn test_retarder_transparent_with_phase_zero(theta1 in 0 as f64..360 as f64,
                                                 theta2 in 0 as f64..360 as f64,) {
        let beam = beam_lin_pol(Angle::Degrees(theta1));
        let retarder = Retarder::new(Angle::Degrees(theta2), Angle::Degrees(0.0));
        let beam_after = retarder.matrix() * beam;
        assert_beam_approx_eq!(beam_after, beam);
    }

    #[test]
    fn test_retarder_transparent_with_phase_2pi(theta1 in 0 as f64..360 as f64,
                                                theta2 in 0 as f64..360 as f64,) {
        let beam = beam_lin_pol(Angle::Degrees(theta1));
        let retarder = Retarder::new(Angle::Degrees(theta2), Angle::Degrees(360.0));
        let beam_after = retarder.matrix() * beam;
        assert_beam_approx_eq!(beam_after, beam);
    }

    #[test]
    fn test_retarder_reduces_to_qwp(theta in 0 as f64..360 as f64) {
        use jones::qwp::QuarterWavePlate;
        let qwp = QuarterWavePlate::new(Angle::Degrees(theta));
        let retarder = Retarder::new(Angle::Degrees(theta), Angle::Degrees(90.0));
        assert_matrix_approx_eq!(qwp.matrix(), retarder.matrix());
    }

    #[test]
    fn test_retarder_reduces_to_hwp(theta in 0 as f64..360 as f64) {
        use jones::hwp::HalfWavePlate;
        let hwp = HalfWavePlate::new(Angle::Degrees(theta));
        let retarder = Retarder::new(Angle::Degrees(theta), Angle::Degrees(180.0));
        assert_matrix_approx_eq!(hwp.matrix(), retarder.matrix());
    }
}
