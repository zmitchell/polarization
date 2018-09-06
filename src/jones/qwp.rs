use na::{Matrix2, Vector2};
use num::complex::Complex;

use super::common::{
    pi, rotate_matrix, Angle, ComplexMatrix, ElementParams, JonesError, JonesMatrix,
    MissingParameter, Result,
};

#[derive(Debug)]
pub struct QuarterWavePlate {
    mat: ComplexMatrix,
}

impl QuarterWavePlate {
    pub fn new(angle: Angle) -> Self {
        let rad = match angle {
            Angle::Radians(ang) => ang,
            Angle::Degrees(ang) => ang.to_radians(),
        };
        let cos_squared = Complex::<f64>::new(rad.cos().powi(2), 0.0);
        let sin_squared = Complex::<f64>::new(rad.sin().powi(2), 0.0);
        let sin_cos = Complex::<f64>::new(rad.sin() * rad.cos(), 0.0);
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

impl From<ElementParams> for Result<QuarterWavePlate> {
    fn from(params: ElementParams) -> Result<QuarterWavePlate> {
        match params.angle {
            Some(angle) => Ok(QuarterWavePlate::new(angle)),
            None => {
                let missing = MissingParameter {
                    typ: "QuarterWavePlate".into(),
                    param: "angle".into(),
                };
                Err(JonesError::MissingParameter(missing))
            }
        }
    }
}

impl JonesMatrix for QuarterWavePlate {
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
    fn test_horizontal_qwp_doesnt_modify_horizontal_beam() {
        let beam = Vector2::new(Complex::<f64>::new(1.0, 0.0), Complex::<f64>::new(0.0, 0.0));
        let qwp = QuarterWavePlate::new(Angle::Degrees(0.0));
        let beam_after = qwp.matrix() * beam;
        assert_beam_approx_eq!(beam, beam_after);
    }

    #[test]
    fn test_vertical_qwp_doesnt_modify_vertical_beam() {
        let beam = Vector2::new(Complex::<f64>::new(0.0, 0.0), Complex::<f64>::new(1.0, 0.0));
        let qwp = QuarterWavePlate::new(Angle::Degrees(90.0));
        let beam_after = qwp.matrix() * beam;
        assert_beam_approx_eq!(beam, beam_after);
    }

    #[test]
    fn test_two_qwps_reflect_polarization() {
        // Horizontal beam, should end up vertical after two QWPs at 45 degrees
        let horiz_beam = Vector2::new(Complex::<f64>::new(1.0, 0.0), Complex::<f64>::new(0.0, 0.0));
        let vert_beam = Vector2::new(Complex::<f64>::new(0.0, 0.0), Complex::<f64>::new(1.0, 0.0));
        let qwp = QuarterWavePlate::new(Angle::Degrees(45.0));
        let beam_after = qwp.matrix() * qwp.matrix() * horiz_beam;
        assert_beam_approx_eq!(vert_beam, beam_after);
    }

    #[test]
    fn test_qwp_circularly_polarizes_beam() {
        // Horizontal beam, should end up right-hand circularly polarized
        let beam = Vector2::new(Complex::<f64>::new(1.0, 0.0), Complex::<f64>::new(0.0, 0.0));
        let qwp = QuarterWavePlate::new(Angle::Degrees(45.0));
        let beam_after = qwp.matrix() * beam;
        // Constant factor out front, e^(i*pi/4) / sqrt(2)
        let prefactor = (Complex::<f64>::i() * pi / 4.0).exp() / (2.0_f64.sqrt());
        // Right-hand circularly polarized beam
        let expected_beam = Vector2::new(prefactor * 1.0, prefactor * (-1.0) * Complex::<f64>::i());
        assert_beam_approx_eq!(beam_after, expected_beam);
    }
}
