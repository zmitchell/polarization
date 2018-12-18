//! Types and definitions used in other modules.
#![macro_use]
use std::error;
pub use std::f64::consts::PI as pi;
use std::fmt;
use std::result;

use na::{Matrix2, Vector2};
use num::complex::Complex;
#[cfg(test)]
use proptest::prelude::*;

/// A more convenient synonym for the type of 2x2 Jones matrices.
pub type ComplexMatrix = Matrix2<Complex<f64>>;

/// A more convenient synonym for the type of 2x1 Jones vectors.
pub type ComplexVector = Vector2<Complex<f64>>;

/// The result type used by `jones`.
///
/// Each error contains an `ErrorKind` to indicate the kind of error encountered.
pub type Result<T> = result::Result<T, JonesError>;

/// The different kinds of errors that may occur inside `polarization`.
#[derive(Debug)]
pub enum JonesError {
    /// An error encountered when the calculated intensity becomes infinite.
    ///
    /// Calculating the intensity involves squaring the components of the Jones vector.
    /// If the components of the vector are large enough, the intensity may become
    /// infinite.
    IntensityTooLarge,

    /// An error encountered when a beam is missing from an optical system.
    NoBeam,

    /// An error encountered when an optical system contains no elements.
    NoElements,

    /// An error not covered by the other kinds of errors.
    Other(String),
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
            JonesError::Other(ref msg) => write!(f, "Other error: {}", msg),
            JonesError::NoBeam => {
                write!(f, "Optical system error: the system does not have a beam")
            }
            JonesError::NoElements => write!(f, "Optical system error: the system has no elements"),
        }
    }
}

/// An angle.
///
/// Angles or phases are more commonly written in radians in physics, but may be more
/// convenient to write in degrees. Furthermore, explicitly denoting the units for
/// angles prevents confusion or mistakes.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Angle {
    Degrees(f64),
    Radians(f64),
}

#[cfg(test)]
impl Arbitrary for Angle {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        prop_oneof![
            (any::<f64>()).prop_map(|x| Angle::Degrees(x)),
            (any::<f64>()).prop_map(|x| Angle::Radians(x)),
        ]
        .boxed()
    }
}

/// The handedness of a circularly polarized beam.
///
/// A circularly polarized beam may either be left- or right-hand circularly polarized,
/// as determined by the right hand rule.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Handedness {
    Left,
    Right,
}

#[cfg(test)]
impl Arbitrary for Handedness {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        prop_oneof![Just(Handedness::Left), Just(Handedness::Right),].boxed()
    }
}

/// The types of polarization handled by Jones calculus.
///
/// While more exotic forms of polarization are possible i.e. vector polarizations, the
/// most common type of polarization are linear, circular, and elliptical. The
/// `Elliptical` variant allows you to specify an arbitrary polarization.
#[derive(Debug, Clone, PartialEq)]
pub enum Polarization {
    Linear(Angle),
    Circular(Handedness),
    Elliptical(Complex<f64>, Complex<f64>),
}

#[cfg(test)]
pub fn any_linear_polarization() -> impl Strategy<Value = Polarization> {
    any::<Angle>().prop_map(|angle| Polarization::Linear(angle))
}

#[cfg(test)]
pub fn any_circular_polarization() -> impl Strategy<Value = Polarization> {
    any::<Handedness>().prop_map(|h| Polarization::Circular(h))
}

#[cfg(test)]
pub fn any_elliptical_polarization() -> impl Strategy<Value = Polarization> {
    (any::<f64>(), any::<f64>(), any::<f64>(), any::<f64>()).prop_map(|(xr, xi, yr, yi)| {
        let x = Complex::new(xr, xi);
        let y = Complex::new(yr, yi);
        Polarization::Elliptical(x, y)
    })
}

#[cfg(test)]
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum PolarizationKind {
    Linear,
    Circular,
    Elliptical,
    Any,
}

#[cfg(test)]
impl Default for PolarizationKind {
    fn default() -> Self {
        PolarizationKind::Any
    }
}

#[cfg(test)]
impl Arbitrary for Polarization {
    type Parameters = PolarizationKind;
    type Strategy = BoxedStrategy<Polarization>;

    fn arbitrary_with(args: Self::Parameters) -> Self::Strategy {
        match args {
            PolarizationKind::Linear => any_linear_polarization().boxed(),
            PolarizationKind::Circular => any_circular_polarization().boxed(),
            PolarizationKind::Elliptical => any_elliptical_polarization().boxed(),
            PolarizationKind::Any => prop_oneof![
                any_linear_polarization(),
                any_circular_polarization(),
                any_elliptical_polarization(),
            ]
            .boxed(),
        }
    }
}

pub trait JonesVector {
    /// Construct a beam from the specified type of polarization. The intensity
    /// of the beam will be 1.
    fn from_polarization(pol: Polarization) -> Self;

    /// The intensity of the beam represented by the Jones vector. Note that
    /// this is simply `V* x V` and is reported as a dimensionless quantity.
    fn intensity(&self) -> Result<f64>;

    /// For the vector `V = (A, B)`, return a new vector with the phase common
    /// to `A` and `B` removed. For the returned vector `V' = (A', B')`, `A'`
    /// will be real, and `B'` will have a phase relative to `A'`.
    fn remove_common_phase(&self) -> Self;

    /// For the vector `V = (A, B)`, remove the phase common to both `A` and `B`
    /// in place. See `remove_common_phase` for more details.
    fn remove_common_phase_mut(&mut self);

    /// Returns the relative phase between the x- and y-components of the vector
    /// in degrees.
    fn relative_phase_deg(&self) -> f64;

    /// Returns the relative phase between the x- and y-components of the vector
    /// in radians.
    fn relative_phase_rad(&self) -> f64;

    /// Returns the 2x1 Jones vector.
    fn vector(&self) -> ComplexVector;

    /// Returns a Jones vector that is the result of passing the current vector through
    /// the provided optical element.
    fn apply_element<T: JonesMatrix>(&self, elem: T) -> Self;

    /// Replace the current Jones vector with the result of passing it through the
    /// provided optical element.
    fn apply_element_mut<T: JonesMatrix>(&mut self, elem: T);

    /// Return the x-component of the Jones vector
    fn x(&self) -> Complex<f64>;

    /// Return the y-component of the Jones vector
    fn y(&self) -> Complex<f64>;
}

/// An ideal coherent light source i.e. an ideal laser beam.
#[derive(Debug, Clone, PartialEq)]
pub struct Beam {
    vec: ComplexVector,
}

impl Beam {
    /// Construct a new beam with arbitrary x- and y-components.
    pub fn new(x: Complex<f64>, y: Complex<f64>) -> Self {
        Beam {
            vec: ComplexVector::new(x, y),
        }
    }

    /// Construct a linearly polarized beam at the angle `angle`, where the angle is
    /// measured counter-clockwise from the x-axis.
    pub fn linear(angle: Angle) -> Self {
        let rad = match angle {
            Angle::Radians(rad) => rad,
            Angle::Degrees(deg) => deg.to_radians(),
        };
        let x = Complex::<f64>::new(rad.cos(), 0.0);
        let y = Complex::<f64>::new(rad.sin(), 0.0);
        Beam {
            vec: ComplexVector::new(x, y),
        }
    }

    /// Construct a circularly polarized beam with the specified handedness.
    pub fn circular(hand: Handedness) -> Self {
        let norm = (1_f64 / 2_f64).sqrt();
        let i = Complex::<f64>::i();
        let x = Complex::<f64>::new(norm, 0.0);
        let y = match hand {
            Handedness::Left => norm * i,
            Handedness::Right => -norm * i,
        };
        Beam {
            vec: ComplexVector::new(x, y),
        }
    }

    /// Construct a beam from an existing 2x1 vector.
    ///
    /// This method is most useful for constructing a `Beam` when you have just received
    /// a `ComplexVector` as the result of some other operation.
    pub fn from_vec(v: ComplexVector) -> Self {
        Beam { vec: v }
    }
}

impl JonesVector for Beam {
    fn from_polarization(pol: Polarization) -> Self {
        use self::Polarization::*;

        match pol {
            Linear(angle) => Beam::linear(angle),
            Circular(hand) => Beam::circular(hand),
            Elliptical(x, y) => Beam::new(x, y),
        }
    }

    fn intensity(&self) -> Result<f64> {
        let conj = ComplexVector::new(self.vec.x.conj(), self.vec.y.conj());
        // The dot product of v* and v should be real, but Rust doesn't know that
        let num = conj.dot(&self.vec).re;
        if num.is_infinite() {
            Err(JonesError::IntensityTooLarge)
        } else {
            Ok(num)
        }
    }

    fn remove_common_phase(&self) -> Self {
        let (x_mag, x_phase) = self.vec.x.to_polar();
        let (y_mag, y_phase) = self.vec.y.to_polar();
        let mut new_y_mag = y_mag;
        let new_x_phase = 0.0;
        let mut new_y_phase = y_phase - x_phase;
        if (y_mag == 0.0) || (y_mag == -0.0) {
            new_y_phase = 0.0;
            new_y_mag = 0.0; // turn -0.0 into 0.0 so that we don't get a phase of pi from the minus sign
        } else {
            // Get the phase in the range [-pi, +pi]
            while new_y_phase < -pi {
                new_y_phase += 2.0 * pi;
            }
            while new_y_phase > pi {
                new_y_phase -= 2.0 * pi;
            }
        }
        let new_y = Complex::from_polar(&new_y_mag, &new_y_phase);
        let new_x = Complex::<f64>::new(x_mag, new_x_phase);
        Beam {
            vec: ComplexVector::new(new_x, new_y),
        }
    }

    fn remove_common_phase_mut(&mut self) {
        let (x_mag, x_phase) = self.vec.x.to_polar();
        let (y_mag, y_phase) = self.vec.y.to_polar();
        let mut new_y_mag = y_mag;
        let new_x_phase = 0.0;
        let mut new_y_phase = y_phase - x_phase;
        if (y_mag == 0.0) || (y_mag == -0.0) {
            new_y_phase = 0.0;
            new_y_mag = 0.0; // turn -0.0 into 0.0 so that we don't get a phase of pi from the minus sign
        } else {
            // Get the phase in the range [-pi, +pi]
            while new_y_phase < -pi {
                new_y_phase += 2.0 * pi;
            }
            while new_y_phase > pi {
                new_y_phase -= 2.0 * pi;
            }
        }
        self.vec.x = Complex::from_polar(&x_mag, &new_x_phase);
        self.vec.y = Complex::from_polar(&new_y_mag, &new_y_phase);
    }

    fn relative_phase_deg(&self) -> f64 {
        let (_, x_phase) = self.vec.x.to_polar();
        let (_, y_phase) = self.vec.y.to_polar();
        y_phase.to_degrees() - x_phase.to_degrees()
    }

    fn relative_phase_rad(&self) -> f64 {
        let (_, x_phase) = self.vec.x.to_polar();
        let (_, y_phase) = self.vec.y.to_polar();
        y_phase - x_phase
    }

    fn vector(&self) -> ComplexVector {
        self.vec
    }

    fn apply_element<T: JonesMatrix>(&self, elem: T) -> Self {
        let vec_after = elem.matrix() * self.vec;
        Beam::from_vec(vec_after)
    }

    fn apply_element_mut<T: JonesMatrix>(&mut self, elem: T) {
        self.vec = elem.matrix() * self.vec;
    }

    fn x(&self) -> Complex<f64> {
        self.vec.x
    }

    fn y(&self) -> Complex<f64> {
        self.vec.y
    }
}

#[cfg(test)]
impl Arbitrary for Beam {
    type Parameters = PolarizationKind;
    type Strategy = BoxedStrategy<Beam>;

    fn arbitrary_with(args: Self::Parameters) -> Self::Strategy {
        match args {
            PolarizationKind::Linear => any_linear_beam().boxed(),
            PolarizationKind::Circular => any_circular_beam().boxed(),
            PolarizationKind::Elliptical => any_elliptical_beam().boxed(),
            PolarizationKind::Any => prop_oneof![
                any_linear_beam(),
                any_circular_beam(),
                any_elliptical_beam(),
            ]
            .boxed(),
        }
    }
}

#[cfg(test)]
pub fn any_linear_beam() -> impl Strategy<Value = Beam> {
    any::<Angle>().prop_map(|angle| Beam::linear(angle))
}

#[cfg(test)]
pub fn any_circular_beam() -> impl Strategy<Value = Beam> {
    any::<Handedness>().prop_map(|h| Beam::circular(h))
}

#[cfg(test)]
pub fn any_elliptical_beam() -> impl Strategy<Value = Beam> {
    (any_complex(), any_complex()).prop_map(|(x, y)| Beam::new(x, y))
}

pub trait JonesMatrix {
    /// Return the optical element rotated by the given angle.
    fn rotated(&self, angle: Angle) -> Self;

    /// Rotate the optical element in-place by the given angle.
    fn rotate(&mut self, angle: Angle);

    /// Returns the 2x2 Jones matrix of the optical element.
    fn matrix(&self) -> ComplexMatrix;
}

/// Rotate an optical element by transforming its Jones matrix.
pub fn rotate_matrix(mat: &ComplexMatrix, angle: &Angle) -> ComplexMatrix {
    let rad = match *angle {
        Angle::Radians(rad) => rad,
        Angle::Degrees(deg) => deg.to_radians(),
    };
    let rot_mat = Matrix2::new(
        Complex::new(rad.cos(), 0_f64),
        Complex::new(rad.sin(), 0_f64),
        Complex::new(-rad.sin(), 0_f64),
        Complex::new(rad.cos(), 0_f64),
    );
    let rot_mat_inv = Matrix2::new(
        Complex::new(rad.cos(), 0_f64),
        Complex::new(-rad.sin(), 0_f64),
        Complex::new(rad.sin(), 0_f64),
        Complex::new(rad.cos(), 0_f64),
    );
    rot_mat_inv * mat * rot_mat
}

#[cfg(test)]
pub fn any_complex() -> impl Strategy<Value = Complex<f64>> {
    (nice_f64(), nice_f64())
        .prop_map(|(x, y)| Complex::new(x, y))
        .boxed()
}

#[cfg(test)]
pub fn nice_f64() -> impl Strategy<Value = f64> {
    (-1e3_f64..1e3_f64)
        .prop_filter("f64 shouldn't be too big", |x| x.abs() < 1e3)
        .prop_filter("f64 shouldn't be too small, but may be zero", |x| {
            (x.abs() > 1e-3) || (*x == 0.0)
        })
}

#[cfg(test)]
pub fn float_angle() -> impl Strategy<Value = f64> {
    (0_f64..180_f64).boxed()
}

#[cfg(test)]
pub fn float_is_well_behaved(x: f64) -> bool {
    let not_nan = !x.is_nan();
    let not_too_small = x > -1e12;
    let not_too_big = x < 1e12;
    return not_nan && not_too_small && not_too_big;
}

#[cfg(test)]
macro_rules! assert_complex_approx_eq {
    ($x:expr, $y:expr) => {
        assert_approx_eq!($x.re, $y.re);
        assert_approx_eq!($x.im, $y.im);
    };
}

#[cfg(test)]
macro_rules! assert_beam_approx_eq {
    ($x:expr, $y:expr) => {
        let vec1 = $x.vector();
        let vec2 = $y.vector();
        assert_complex_approx_eq!(vec1[0], vec2[0]);
        assert_complex_approx_eq!(vec1[1], vec2[1]);
    };
}

#[cfg(test)]
macro_rules! assert_matrix_approx_eq {
    ($x:expr, $y:expr) => {
        assert_complex_approx_eq!($x[(0, 0)], $y[(0, 0)]);
        assert_complex_approx_eq!($x[(0, 1)], $y[(0, 1)]);
        assert_complex_approx_eq!($x[(1, 0)], $y[(1, 0)]);
        assert_complex_approx_eq!($x[(1, 1)], $y[(1, 1)]);
    };
}

#[cfg(test)]
macro_rules! prop_assert_approx_eq {
    ($x:expr, $y:expr) => {
        prop_assert!(($x - $y).abs() < 1e-6)
    };
}

#[cfg(test)]
macro_rules! prop_assert_complex_approx_eq {
    ($x:expr, $y:expr) => {
        prop_assert_approx_eq!($x.re, $y.re);
        prop_assert_approx_eq!($x.im, $y.im);
    };
}

#[cfg(test)]
macro_rules! prop_assert_beam_approx_eq {
    ($x:expr, $y:expr) => {
        let vec1 = $x.vector();
        let vec2 = $y.vector();
        prop_assert_complex_approx_eq!(vec1[0], vec2[0]);
        prop_assert_complex_approx_eq!(vec1[1], vec2[1]);
    };
}

#[cfg(test)]
macro_rules! prop_assert_matrix_approx_eq {
    ($x:expr, $y:expr) => {
        prop_assert_complex_approx_eq!($x[(0, 0)], $y[(0, 0)]);
        prop_assert_complex_approx_eq!($x[(0, 1)], $y[(0, 1)]);
        prop_assert_complex_approx_eq!($x[(1, 0)], $y[(1, 0)]);
        prop_assert_complex_approx_eq!($x[(1, 1)], $y[(1, 1)]);
    };
}

// Beam tests
#[cfg(test)]
proptest! {
   #[test]
   fn test_intensity_one_real_component(n: f64) {
       let beam = Beam::new(
           Complex::new(n, 0_f64),
           Complex::new(0_f64, 0_f64),
       );
       let intensity = beam.intensity();
       prop_assume!(intensity.is_ok());
       prop_assert_approx_eq!(n.powi(2), intensity.unwrap());
   }

   #[test]
   fn test_intensity_is_never_negative(beam: Beam) {
       prop_assume!(beam.intensity().is_ok());
       let intensity = beam.intensity().unwrap();
       prop_assert!(intensity >= 0.0);
   }

   #[test]
   fn test_intensity_mag_with_complex(x in any_complex(), y in any_complex()) {
       let xr = x.re;
       let xi = x.im;
       let yr = y.re;
       let yi = y.im;
       let by_hand = xr.powi(2) + xi.powi(2) + yr.powi(2) + yi.powi(2);
       let beam = Beam::new(x, y);
       prop_assume!(beam.intensity().is_ok());
       prop_assert_approx_eq!(beam.intensity().unwrap(), by_hand);
   }

   #[test]
   fn test_common_phase_preserves_x_mag(beam: Beam) {
       let old_x_mag = beam.x().norm();
       let new_beam = beam.remove_common_phase();
       let new_x_mag = new_beam.x().norm();
       prop_assert_eq!(old_x_mag, new_x_mag);
   }

   #[test]
   fn test_common_phase_preserves_y_mag(beam: Beam) {
       let old_y_mag = beam.y().norm();
       let new_beam = beam.remove_common_phase();
       let new_y_mag = new_beam.y().norm();
       prop_assert!((old_y_mag - new_y_mag).abs() < 1e-6);
   }

   #[test]
   fn test_common_phase_zeroes_x_phase(beam: Beam) {
       let new_beam = beam.remove_common_phase();
       let new_phase = new_beam.x().arg();
       prop_assert!(new_phase.abs() < 1e-6);
   }

   #[test]
   fn test_common_phase_correct_y_phase(beam: Beam) {
       let old_y_phase = beam.y().arg();
       let old_x_phase = beam.x().arg();
       let new_beam = beam.remove_common_phase();
       let new_y_phase = new_beam.y().arg();
       let mut expected_y_phase = old_y_phase - old_x_phase;
       if beam.y().norm() == 0.0 {
           // Regardless of the difference between old_y_phase and old_x_phase, if the magnitude of
           // y is +0.0, the new phase will be wiped out. The complex number is stored as real and
           // imaginary parts computed as mag*cos(phase) and mag*sin(phase) respectively, so when
           // the magnitude is zero, mag*sin(phase) is also zero. This makes it such that the phase
           // you use to construct the complex number will not be the phase returned by
           // beam.y().arg().
           expected_y_phase = 0.0;
       } else {
           // Get the phase in the range [-pi, pi]
           while expected_y_phase < -pi {
               expected_y_phase = expected_y_phase + 2.0 * pi;
           }
           while expected_y_phase > pi {
               expected_y_phase = expected_y_phase - 2.0 * pi;
           }
       }
       // If the magnitude is zero, the phase will also be set to zero due to how
       // Complex::from_polar is implemented, so the phases won't match up. In this
       // test I'll check that the magnitude is zero if the phases don't match.
       if (expected_y_phase - new_y_phase).abs() > 1e-6 {
           let new_norm = new_beam.y().norm();
           if new_norm.abs() == 0.0 {
               prop_assert!(new_y_phase.abs() < 1e-6);
           }
       } else {
           prop_assert!((expected_y_phase - new_y_phase).abs() < 1e-6);
       }
   }

   #[test]
   fn test_common_phase_preserves_intensity(beam: Beam) {
       let old_intensity = beam.intensity();
       prop_assume!(old_intensity.is_ok());
       let new_beam = beam.remove_common_phase();
       let new_intensity = new_beam.intensity();
       prop_assert!(new_intensity.is_ok());
       prop_assert!((new_intensity.unwrap() - old_intensity.unwrap()).abs() < 1e-6);
   }

    #[test]
    fn test_common_phase_mut_preserves_x_mag(beam: Beam) {
        let old_x_mag = beam.x().norm();
        beam.remove_common_phase();
        let new_x_mag = beam.x().norm();
        prop_assert_eq!(old_x_mag, new_x_mag);
    }

    #[test]
    fn test_common_phase_mut_preserves_y_mag(beam: Beam) {
        let old_y_mag = beam.y().norm();
        beam.remove_common_phase();
        let new_y_mag = beam.y().norm();
        prop_assert!((old_y_mag - new_y_mag).abs() < 1e-6);
    }

    #[test]
    fn test_common_phase_mut_zeroes_x_phase(mut beam: Beam) {
        beam.remove_common_phase_mut();
        let new_phase = beam.x().arg();
        prop_assert!(new_phase.abs() < 1e-6);
    }

    #[test]
    fn test_common_phase_mut_correct_y_phase(mut beam: Beam) {
        let old_y_phase = beam.y().arg();
        let old_x_phase = beam.x().arg();
        beam.remove_common_phase_mut();
        let new_y_phase = beam.y().arg();
        let mut expected_y_phase = old_y_phase - old_x_phase;
        if beam.y().norm() == 0.0 {
           // Regardless of the difference between old_y_phase and old_x_phase, if the magnitude of
           // y is +0.0, the new phase will be wiped out. The complex number is stored as real and
           // imaginary parts computed as mag*cos(phase) and mag*sin(phase) respectively, so when
           // the magnitude is zero, mag*sin(phase) is also zero. This makes it such that the phase
           // you use to construct the complex number will not be the phase returned by
           // beam.y().arg().
            expected_y_phase = 0.0;
        } else {
            // Get the phase in the range [-pi, pi]
            while expected_y_phase < - pi {
                expected_y_phase = expected_y_phase + 2.0 * pi;
            }
            while expected_y_phase > pi {
                expected_y_phase = expected_y_phase - 2.0 * pi;
            }
        }
        // If the magnitude is zero, the phase will also be set to zero due to how
        // Complex::from_polar is implemented, so the phases won't match up. In this
        // test I'll check that the magnitude is zero if the phases don't match.
        if (expected_y_phase - new_y_phase).abs() > 1e-6 {
            let new_norm = beam.y().norm();
            if new_norm.abs() == 0.0 {
                prop_assert!(new_y_phase.abs() < 1e-6);
            }
        } else {
            prop_assert!((expected_y_phase - new_y_phase).abs() < 1e-6);
        }
    }

    #[test]
    fn test_common_phase_mut_preserves_intensity(beam: Beam) {
        let old_intensity = beam.intensity();
        prop_assume!(old_intensity.is_ok());
        beam.remove_common_phase();
        let new_intensity = beam.intensity();
        prop_assert!(new_intensity.is_ok());
        prop_assert!((new_intensity.unwrap() - old_intensity.unwrap()).abs() < 1e-6);
    }


   #[test]
   fn test_relative_phase(x in any_complex(), y in any_complex()) {
       let expected_rad = y.arg() - x.arg();
       let expected_deg = expected_rad.to_degrees();
       let beam = Beam::new(x, y);
       prop_assert_approx_eq!(expected_rad, beam.relative_phase_rad());
       prop_assert_approx_eq!(expected_deg, beam.relative_phase_deg());
   }

    #[test]
    fn test_construct_beam_from_xy(x in any_complex(), y in any_complex()) {
        let beam = Beam::new(x, y);
        prop_assert_eq!(beam.x(), x);
        prop_assert_eq!(beam.y(), y);
    }

    #[test]
    fn test_rotate_360_degrees_returns_original(m00 in any_complex(),
                                                m01 in any_complex(),
                                                m10 in any_complex(),
                                                m11 in any_complex()) {
        let mat = Matrix2::new(m00, m01, m10, m11);
        let angle = Angle::Degrees(360_f64);
        let rotated = rotate_matrix(&mat, &angle);
        prop_assert_matrix_approx_eq!(mat, rotated);
    }

   #[test]
   fn test_rotate_2pi_rad_returns_original(m00 in any_complex(),
                                           m01 in any_complex(),
                                           m10 in any_complex(),
                                           m11 in any_complex()) {
       let mat = Matrix2::new(m00, m01, m10, m11);
       let angle = Angle::Radians(2.0 * pi);
       let rotated = rotate_matrix(&mat, &angle);
       prop_assert_matrix_approx_eq!(mat, rotated);
   }

   #[test]
   fn test_rotate_0_degrees_returns_original(m00 in any_complex(),
                                             m01 in any_complex(),
                                             m10 in any_complex(),
                                             m11 in any_complex()) {
       let mat = Matrix2::new(m00, m01, m10, m11);
       let angle = Angle::Degrees(0_f64);
       let rotated = rotate_matrix(&mat, &angle);
       prop_assert_matrix_approx_eq!(mat, rotated);
   }

   #[test]
   fn test_rotate_0_rad_returns_original(m00 in any_complex(),
                                         m01 in any_complex(),
                                         m10 in any_complex(),
                                         m11 in any_complex()) {
       let mat = Matrix2::new(m00, m01, m10, m11);
       let angle = Angle::Radians(0_f64);
       let rotated = rotate_matrix(&mat, &angle);
       prop_assert_matrix_approx_eq!(mat, rotated);
   }

}
