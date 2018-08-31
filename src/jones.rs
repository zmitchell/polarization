use std::error;
use std::f64::consts::PI as pi;
use std::fmt;
use std::result;

use na::{Matrix2, Vector2};
use num::complex::Complex;
use proptest;
use proptest::num::f64::{NEGATIVE, POSITIVE, ZERO};
use proptest::prelude::*;

use core::{Angle, BeamPol, Handedness};

type ComplexMatrix = Matrix2<Complex<f64>>;
type Result<T> = result::Result<T, JonesError>;

#[derive(Debug)]
pub enum JonesError {
    IntensityTooLarge,
    IntensityTooSmall,
    MissingParameter(MissingParameter),
    Other(String),
}

#[derive(Debug)]
pub struct MissingParameter {
    pub typ: String,
    pub param: String,
}

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
            JonesError::MissingParameter(ref missing_param) => write!(
                f,
                "Missing parameter: {} requires parameter '{}'",
                missing_param.typ, missing_param.param
            ),
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
        let conj = Vector2::new(self.x.conj(), self.y.conj());
        // The dot product of v* and v should be real, but Rust doesn't know that
        let num = conj.dot(self).re;
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
        Vector2::new(Complex::from_polar(&x_mag, &0.0_f64), y_rel)
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
        }
        Angle::Degrees(deg) => {
            x = deg.to_radians().cos();
            y = deg.to_radians().sin();
        }
    }
    Vector2::new(Complex::new(x, 0_f64), Complex::new(y, 0_f64))
}

pub trait JonesMatrix {
    /// Produce the Jones matrix after rotating the element by the given angle.
    fn rotated(&self, angle: Angle) -> ComplexMatrix;

    /// Rotate the element in-place by the given angle.
    fn rotate(&mut self, angle: Angle);

    /// Returns the Jones matrix of the optical element.
    fn matrix(&self) -> ComplexMatrix;
}

/// Returns the matrix of an optical element after it has been rotated around the optical axis by
/// the given angle.
pub fn rotate_matrix(mat: &ComplexMatrix, angle: &Angle) -> ComplexMatrix {
    let rad = match angle {
        &Angle::Radians(rad) => rad,
        &Angle::Degrees(deg) => deg.to_radians(),
    };
    let rot_mat = Matrix2::new(
        Complex::new(rad.cos(), 0 as f64),
        Complex::new(rad.sin(), 0 as f64),
        Complex::new(-rad.sin(), 0 as f64),
        Complex::new(rad.cos(), 0 as f64),
    );
    let rot_mat_inv = Matrix2::new(
        Complex::new(rad.cos(), 0 as f64),
        Complex::new(-rad.sin(), 0 as f64),
        Complex::new(rad.sin(), 0 as f64),
        Complex::new(rad.cos(), 0 as f64),
    );
    rot_mat_inv * mat * rot_mat
}

#[derive(Debug)]
pub enum OpticalElement {
    Polarizer,
    QuarterWavePlate,
    HalfWavePlate,
    PolarizationRotator,
    Retarder,
    DieletricReflection,
    MetalReflection,
}

#[derive(Debug)]
pub struct ElementParams {
    pub angle: Option<Angle>,
    pub incident_angle: Option<Angle>,
    pub azimuthal_angle: Option<Angle>,
    pub refractive_index: Option<f64>,
    pub extinction_coefficient: Option<f64>,
}

#[derive(Debug)]
pub struct Polarizer {
    angle: Angle,
    mat: ComplexMatrix,
}

impl Polarizer {
    pub fn new(angle: Angle) -> Self {
        let x_pol = Matrix2::new(
            Complex::new(1 as f64, 0 as f64),
            Complex::new(0 as f64, 0 as f64),
            Complex::new(0 as f64, 0 as f64),
            Complex::new(0 as f64, 0 as f64),
        );
        let mat_at_angle = rotate_matrix(&x_pol, &angle);
        Polarizer {
            angle: angle,
            mat: mat_at_angle,
        }
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
        self.mat.clone()
    }
}

#[derive(Debug)]
pub struct QuarterWavePlate {
    angle: Angle,
    mat: ComplexMatrix,
}

#[derive(Debug)]
pub struct HalfWavePlate {
    angle: Angle,
    mat: ComplexMatrix,
}

#[derive(Debug)]
pub struct Retarder {
    angle: Angle,
    phase: Angle,
    mat: ComplexMatrix,
}

#[derive(Debug)]
pub struct PolarizationRotator {
    angle: Angle,
    mat: ComplexMatrix,
}

#[derive(Debug)]
pub struct DielectricReflection {
    refractive_index: f64,
    incident_angle: Angle,
    azimuthal_angle: Angle,
    mat: ComplexMatrix,
}

#[derive(Debug)]
pub struct MetalReflection {
    incident_angle: Angle,
    azimuthal_angle: Angle,
    mat: ComplexMatrix,
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

macro_rules! assert_complex_approx_eq {
    ($x:expr, $y:expr) => {
        assert_approx_eq!($x.re, $y.re);
        assert_approx_eq!($x.im, $y.im);
    };
}

macro_rules! assert_matrix_approx_eq {
    ($x:expr, $y:expr) => {
        assert_complex_approx_eq!($x[(0, 0)], $y[(0, 0)]);
        assert_complex_approx_eq!($x[(0, 1)], $y[(0, 1)]);
        assert_complex_approx_eq!($x[(1, 0)], $y[(1, 0)]);
        assert_complex_approx_eq!($x[(1, 1)], $y[(1, 1)]);
    };
}

proptest!{

    #[test]
    fn test_intensity_one_real_component(n in well_behaved_doubles()) {
        let beam = Vector2::new(
            Complex::new(n, 0_f64),
            Complex::new(0_f64, 0_f64),
        );
        let intensity = beam.intensity().unwrap();
        assert_approx_eq!(n.powi(2), intensity);
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

    #[test]
    fn test_rotate_360_degrees_returns_original(m00 in well_behaved_complexes(),
                                                m01 in well_behaved_complexes(),
                                                m10 in well_behaved_complexes(),
                                                m11 in well_behaved_complexes()) {
        let mat = Matrix2::new(m00, m01, m10, m11);
        let angle = Angle::Degrees(360 as f64);
        let rotated = rotate_matrix(&mat, &angle);
        assert_matrix_approx_eq!(mat, rotated);
    }

    #[test]
    fn test_rotate_2pi_rad_returns_original(m00 in well_behaved_complexes(),
                                            m01 in well_behaved_complexes(),
                                            m10 in well_behaved_complexes(),
                                            m11 in well_behaved_complexes()) {
        let mat = Matrix2::new(m00, m01, m10, m11);
        let angle = Angle::Radians(2.0 * pi);
        let rotated = rotate_matrix(&mat, &angle);
        assert_matrix_approx_eq!(mat, rotated);
    }

    #[test]
    fn test_rotate_0_degrees_returns_original(m00 in well_behaved_complexes(),
                                              m01 in well_behaved_complexes(),
                                              m10 in well_behaved_complexes(),
                                              m11 in well_behaved_complexes()) {
        let mat = Matrix2::new(m00, m01, m10, m11);
        let angle = Angle::Degrees(0 as f64);
        let rotated = rotate_matrix(&mat, &angle);
        assert_matrix_approx_eq!(mat, rotated);
    }

    #[test]
    fn test_rotate_0_rad_returns_original(m00 in well_behaved_complexes(),
                                          m01 in well_behaved_complexes(),
                                          m10 in well_behaved_complexes(),
                                          m11 in well_behaved_complexes()) {
        let mat = Matrix2::new(m00, m01, m10, m11);
        let angle = Angle::Radians(0 as f64);
        let rotated = rotate_matrix(&mat, &angle);
        assert_matrix_approx_eq!(mat, rotated);
    }

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

#[cfg(test)]
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

#[cfg(test)]
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
