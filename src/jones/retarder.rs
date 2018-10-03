//! An optical retarder.
//!
//! A [waveplate](https://en.wikipedia.org/wiki/Waveplate) is an optical device that is used to
//! introduce a phase difference between the two perpendicular components of a beam's
//! polarization. Common waveplates are half-wave plates and quarter-wave plates, which introduce
//! phase differences of `pi` and `pi/2` respectively. A waveplate that introduces an arbitrary
//! phase difference is commonly referred to as a retarder. The phase difference introduced by the
//! retarder is commonly quoted as an angle i.e. 90 degrees or `pi` radians.
use na::Matrix2;
use num::complex::Complex;

use super::common::{rotate_matrix, Angle, ComplexMatrix, JonesMatrix};

/// An ideal optical retarder that can introduce an arbitrary phase.
#[derive(Debug, Copy, Clone)]
pub struct Retarder {
    mat: ComplexMatrix,
}

impl Retarder {
    /// Constructs a new retarder.
    ///
    /// The constructor takes two arguments, both of which are angles. The first angle describes
    /// the orientation of the fast-axis of the retarder. The second angle describes the phase
    /// difference that the retarder introduces.
    pub fn new(angle: Angle, phase: Angle) -> Self {
        let angle_rad = match angle {
            Angle::Degrees(deg) => deg.to_radians(),
            Angle::Radians(rad) => rad,
        };
        let phase_rad = match phase {
            Angle::Degrees(deg) => deg.to_radians(),
            Angle::Radians(rad) => rad,
        };
        let sin = Complex::new(angle_rad.sin(), 0_f64);
        let cos = Complex::new(angle_rad.cos(), 0_f64);
        let sin_2 = Complex::new(angle_rad.sin().powi(2), 0_f64);
        let cos_2 = Complex::new(angle_rad.cos().powi(2), 0_f64);
        let xi = (Complex::<f64>::i() * phase_rad).exp();
        let a = cos_2 + xi * sin_2;
        let b = sin * cos - xi * sin * cos;
        let c = sin * cos - xi * sin * cos;
        let d = sin_2 + xi * cos_2;
        let mat = Matrix2::new(a, b, c, d);
        Retarder { mat }
    }
}

impl JonesMatrix for Retarder {
    /// Returns the element rotated counter-clockwise by `angle`.
    fn rotated(&self, angle: Angle) -> Self {
        Retarder {
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
mod tests {
    use super::*;
    use jones::common::{Beam, JonesVector};

    proptest! {
        #[test]
        fn test_retarder_transparent_with_phase_zero(theta1 in 0_f64..360_f64,
                                                     theta2 in 0_f64..360_f64,) {
            let beam = Beam::linear(Angle::Degrees(theta1));
            let retarder = Retarder::new(Angle::Degrees(theta2), Angle::Degrees(0.0));
            let beam_after = beam.apply_element(retarder);
            assert_beam_approx_eq!(beam_after, beam);
        }

        #[test]
        fn test_retarder_transparent_with_phase_2pi(theta1 in 0_f64..360_f64,
                                                    theta2 in 0_f64..360_f64,) {
            let beam = Beam::linear(Angle::Degrees(theta1));
            let retarder = Retarder::new(Angle::Degrees(theta2), Angle::Degrees(360.0));
            let beam_after = beam.apply_element(retarder);
            assert_beam_approx_eq!(beam_after, beam);
        }

        #[test]
        fn test_retarder_reduces_to_qwp(theta in 0_f64..360_f64) {
            use jones::qwp::QuarterWavePlate;
            let qwp = QuarterWavePlate::new(Angle::Degrees(theta));
            let retarder = Retarder::new(Angle::Degrees(theta), Angle::Degrees(90.0));
            assert_matrix_approx_eq!(qwp.matrix(), retarder.matrix());
        }

        #[test]
        fn test_retarder_reduces_to_hwp(theta in 0_f64..360_f64) {
            use jones::hwp::HalfWavePlate;
            let hwp = HalfWavePlate::new(Angle::Degrees(theta));
            let retarder = Retarder::new(Angle::Degrees(theta), Angle::Degrees(180.0));
            assert_matrix_approx_eq!(hwp.matrix(), retarder.matrix());
        }
    }
}
