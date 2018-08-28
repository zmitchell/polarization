use std::error;
use std::f64::consts::PI as pi;
use std::fmt;
use std::result;

use na::{Matrix2, Vector2};
use num::complex::Complex;
use proptest;
use proptest::prelude::*;
use proptest::num::f64::{POSITIVE, NEGATIVE, ZERO};

use core::{Angle, BeamPol, Handedness};

#[derive(Debug)]
pub enum JonesError {
    IntensityTooLarge,
    IntensityTooSmall,
    Other(String),
}

type Result<T> = result::Result<T, JonesError>;

impl error::Error for JonesError {
    fn description(&self) -> &str {
        ""
    }

    fn cause(&self) -> Option<&error::Error> {
        None
    }
}

impl fmt::Display for JonesError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            JonesError::IntensityTooLarge => write!(f, "Intensity error: Intensity is too large"),
            JonesError::IntensityTooSmall => write!(f, "Intensity error: Intensity is too small"),
            JonesError::Other(ref msg) => write!(f, "Other error: {}", msg),
        }
    }
}


pub trait JonesVector {
    /// The intensity of the beam represented by the Jones vector. Note that
    /// this is simply `V* x V`.
    fn intensity(&self) -> Result<f64>;

    /// For the vector `V = (A, B)`, return a new vector with the phase common
    /// to `A` and `B` removed. For the returned vector `V' = (A', B')`, `A'`
    /// will be real, and `B'` will have some phase relative to `A'`.
    fn remove_common_phase(&self) -> Self;

    /// For the vector `V = (A, B)`, remove the phase common to both `A` and `B`
    /// in place. See `remove_common_phase`.
    fn remove_common_phase_mut(&mut self);

    /// Returns the relative phase between the x- and y-components of the vector
    /// in degrees.
    fn relative_phase_deg(&self) -> f64;

    /// Returns the relative phase between the x- and y-components of the vector
    /// in radians.
    fn relative_phase_rad(&self) -> f64;
}

impl JonesVector for Vector2<Complex<f64>> {
    fn intensity(&self) -> Result<f64> {
        let conj = Vector2::new(
            self.x.conj(),
            self.y.conj(),
        );
        let num = conj.dot(self).norm_sqr().sqrt().sqrt();
        if num.is_infinite() {
            Err(JonesError::IntensityTooLarge)
        } else if (num.abs() < 1e-12) && (num.abs() > 0.0) {
            Err(JonesError::IntensityTooSmall)
        } else {
            Ok(num)
        }
    }

    fn remove_common_phase(&self) -> Self {
        let (x_mag, x_phase) = self.x.to_polar();
        let (y_mag, y_phase) = self.y.to_polar();
        let rel_phase = (y_phase - x_phase) % (2 as f64 * pi);
        let y_rel = Complex::from_polar(&y_mag, &rel_phase);
        Vector2::new(
            Complex::from_polar(&x_mag, &0.0_f64),
            y_rel,
        )
    }

    fn remove_common_phase_mut(&mut self) {
        let (x_mag, x_phase) = self.x.to_polar();
        let (y_mag, y_phase) = self.y.to_polar();
        let rel_phase = (y_phase - x_phase) % (2 as f64 * pi);
        self.x = Complex::from_polar(&x_mag, &0.0_f64);
        self.y = Complex::from_polar(&y_mag, &rel_phase);
    }

    fn relative_phase_deg(&self) -> f64 {
        let (_, x_phase) = self.x.to_polar();
        let (_, y_phase) = self.y.to_polar();
        y_phase.to_degrees() - x_phase.to_degrees()
    }

    fn relative_phase_rad(&self) -> f64 {
        let (_, x_phase) = self.x.to_polar();
        let (_, y_phase) = self.y.to_polar();
        y_phase - x_phase
    }
}

/// Produces the Jones vector for a linearly polarized beam that
/// is parallel to the x-axis.
pub fn beam_horizontal() -> Vector2<Complex<f64>> {
    // v = (1, 0)
    Vector2::new(
        Complex::new(1.0_f64, 0.0_f64),
        Complex::new(0.0_f64, 0.0_f64),
    )
}

/// Produces the Jones vector for a linearly polarized beam that
/// is perpendicular to the x-axis.
pub fn beam_vertical() -> Vector2<Complex<f64>> {
    // v = (0, 1)
    Vector2::new(
        Complex::new(0.0_f64, 0.0_f64),
        Complex::new(1.0_f64, 0.0_f64),
    )
}

/// Produces the Jones vector for a beam with left-handed circular polarization.
pub fn beam_left_circular() -> Vector2<Complex<f64>> {
    // v = (1/sqrt(2)) * (1, i)
    let x = Complex::new(1.0_f64 / 2.0_f64.sqrt(), 0.0_f64);
    let y = Complex::i() / 2.0_f64.sqrt();
    Vector2::new(x, y)
}

/// Produces the Jones vector for a beam with right-handed circular polarization.
pub fn beam_right_circular() -> Vector2<Complex<f64>> {
    // v = (1/sqrt(2)) * (1, i)
    let x = Complex::new(1.0_f64 / 2.0_f64.sqrt(), 0.0_f64);
    let y = -Complex::i() / 2.0_f64.sqrt();
    Vector2::new(x, y)
}

/// Produces a linearly polarized beam at the given angle.
pub fn beam_lin_pol(angle: Angle) -> Vector2<Complex<f64>> {
    let x: f64;
    let y: f64;
    match angle {
        Angle::Radians(rad) => {
            x = rad.cos();
            y = rad.sin();
        },
        Angle::Degrees(deg) => {
            x = deg.to_radians().cos();
            y = deg.to_radians().sin();
        },
    }
    Vector2::new(
        Complex::new(x, 0_f64),
        Complex::new(y, 0_f64),
    )
}

prop_compose! {
    fn well_behaved_doubles()(x in (POSITIVE | NEGATIVE | ZERO).prop_filter(
            "Floats should be zero, or between 1e-12 and 1e12",
            |f| (f == &0_f64) || ((f.abs() < 1e12_f64) && (f.abs() > 1e-12_f64)))
    ) -> f64 {
        x
    }
}

prop_compose! {
    fn well_behaved_complexes()(x in well_behaved_doubles(),
                                y in well_behaved_doubles(),
    ) -> Complex<f64> {
        Complex::new(x, y)
    }
}

// prop_compose! {
//     fn well_behaved_jones_vector()(x in well_behaved_complexes(),
//                                    y in well_behaved_complexes(),
//     ) -> Vector2<Complex<f64>> {
//         Vector2::new(0, 0)
//     }
// }

proptest!{

    #[test]
    fn test_intensity_one_real_component(n in well_behaved_doubles()) {
        let beam = Vector2::new(
            Complex::new(n, 0_f64),
            Complex::new(0_f64, 0_f64),
        );
        let intensity_attempt = beam.intensity();
        assert!(intensity_attempt.is_ok());
        let intensity = intensity_attempt.unwrap();
        assert_approx_eq!(n, intensity);
    }

    #[test]
    fn test_intensity_is_never_negative(x in well_behaved_complexes(),
                                        y in well_behaved_complexes()) {
        let v = Vector2::new(x, y);
        prop_assume!(v.intensity().is_ok());
        let intensity = v.intensity().unwrap();
        assert!(intensity >= 0.0);
    }

    #[test]
    fn test_intensity_mag_with_complex(x in well_behaved_complexes(),
                                       y in well_behaved_complexes()) {
        let xr = x.re;
        let xi = x.im;
        let yr = y.re;
        let yi = y.im;
        let by_hand = xr.powi(2) + xi.powi(2) + yr.powi(2) + yi.powi(2);
        let v = Vector2::new(x, y);
        prop_assume!(v.intensity().is_ok());
        assert_approx_eq!(v.intensity().unwrap(), by_hand);
    }

    #[test]
    fn test_common_phase(x in well_behaved_complexes(),
                         y in well_behaved_complexes()) {
        let y_phase_new = y.arg() - x.arg();
        let beam = Vector2::new(x, y);
        let new_beam = beam.remove_common_phase();
        assert_approx_eq!(new_beam.x.arg(), 0.0 as f64);
        assert_approx_eq!(new_beam.y.arg(), y_phase_new % (2_f64 * pi));
        //
    }

    #[test]
    fn test_common_phase_mut(x in well_behaved_complexes(),
                             y in well_behaved_complexes()) {
        let y_phase_new = y.arg() - x.arg();
        let mut beam = Vector2::new(x, y);
        beam.remove_common_phase_mut();
        assert_approx_eq!(beam.x.arg(), 0.0 as f64);
        assert_approx_eq!(beam.y.arg(), y_phase_new % (2_f64 * pi));
    }

    #[test]
    fn test_relative_phase(x in well_behaved_complexes(),
                           y in well_behaved_complexes()) {
        let expected_rad = y.arg() - x.arg();
        let expected_deg = expected_rad.to_degrees();
        let beam = Vector2::new(x, y);
        assert_approx_eq!(expected_rad, beam.relative_phase_rad());
        assert_approx_eq!(expected_deg, beam.relative_phase_deg());
    }
}
