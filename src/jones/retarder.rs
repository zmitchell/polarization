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
#[cfg(test)]
use proptest::prelude::*;

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

// The `Angle` supplied to each variant is the value of the other parameter i.e.
// `RetarderParam::Phase(angle)` would fix the angle of the retarder at `angle` while supplying
// arbitrary phases.
#[cfg(test)]
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum RetarderParam {
    Phase(Angle),
    Angle(Angle),
    All,
}

#[cfg(test)]
impl Default for RetarderParam {
    fn default() -> Self {
        RetarderParam::All
    }
}

#[cfg(test)]
impl Arbitrary for Retarder {
    type Parameters = RetarderParam;
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(args: Self::Parameters) -> Self::Strategy {
        match args {
            RetarderParam::Phase(angle) => any_retarder_fixed_angle(angle).boxed(),
            RetarderParam::Angle(phase) => any_retarder_fixed_phase(phase).boxed(),
            RetarderParam::All => any_retarder().boxed(),
        }
        .boxed()
    }
}

#[cfg(test)]
fn any_retarder_fixed_angle(angle: Angle) -> impl Strategy<Value = Retarder> {
    any::<Angle>().prop_map(move |phase| Retarder::new(angle, phase))
}

#[cfg(test)]
fn any_retarder_fixed_phase(phase: Angle) -> impl Strategy<Value = Retarder> {
    any::<Angle>().prop_map(move |angle| Retarder::new(angle, phase))
}

#[cfg(test)]
fn any_retarder() -> impl Strategy<Value = Retarder> {
    (any::<Angle>(), any::<Angle>()).prop_map(|(angle, phase)| Retarder::new(angle, phase))
}

#[cfg(test)]
mod tests {
    use super::*;
    use jones::common::{Beam, JonesVector};

    proptest! {
        #[test]
        fn test_retarder_transparent_with_phase_zero(beam: Beam, angle: Angle) {
            let ret = Retarder::new(angle, Angle::Degrees(0.0));
            let beam_after = beam.apply_element(ret);
            prop_assert_beam_approx_eq!(beam_after, beam);
        }

        #[test]
        fn test_retarder_transparent_with_phase_2pi(beam: Beam, angle: Angle) {
            let ret = Retarder::new(angle, Angle::Degrees(360.0));
            let beam_after = beam.apply_element(ret);
            prop_assert_beam_approx_eq!(beam_after, beam);
        }

        #[test]
        fn test_retarder_reduces_to_qwp(angle: Angle) {
            use jones::qwp::QuarterWavePlate;
            let qwp = QuarterWavePlate::new(angle);
            prop_assume!(!qwp.matrix()[(0,0)].is_nan());
            prop_assume!(!qwp.matrix()[(0,1)].is_nan());
            prop_assume!(!qwp.matrix()[(1,0)].is_nan());
            prop_assume!(!qwp.matrix()[(1,1)].is_nan());
            let retarder = Retarder::new(angle, Angle::Degrees(90.0));
            prop_assert_matrix_approx_eq!(qwp.matrix(), retarder.matrix());
        }

        #[test]
        fn test_retarder_reduces_to_hwp(angle: Angle) {
            use jones::hwp::HalfWavePlate;
            let hwp = HalfWavePlate::new(angle);
            prop_assume!(!hwp.matrix()[(0,0)].is_nan());
            prop_assume!(!hwp.matrix()[(0,1)].is_nan());
            prop_assume!(!hwp.matrix()[(1,0)].is_nan());
            prop_assume!(!hwp.matrix()[(1,1)].is_nan());
            let retarder = Retarder::new(angle, Angle::Degrees(180.0));
            prop_assert_matrix_approx_eq!(hwp.matrix(), retarder.matrix());
        }

        #[test]
        fn test_retarder_preserves_intensity(beam: Beam, ret: Retarder) {
            let intensity_before = beam.intensity();
            prop_assume!(intensity_before.is_ok());
            let intensity_after = beam.apply_element(ret).intensity();
            prop_assume!(intensity_after.is_ok());
            prop_assert_approx_eq!(intensity_after.unwrap(), intensity_before.unwrap());
        }
    }
}
