//! Types and definitions used in other modules.
#![macro_use]
use std::error;
pub use std::f64::consts::PI as pi;
use std::fmt;
use std::result;

use na::{Matrix2, Vector2};
use num::complex::Complex;
#[cfg(test)]
use proptest::num::f64::{NEGATIVE, POSITIVE, ZERO};
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

/// The types of polarization handled by Jones calculus.
///
/// While more exotic forms of polarization are possible i.e. vector polarizations, the
/// most common type of polarization are linear, circular, and elliptical. The
/// `Elliptical` variant allows you to specify an arbitrary polarization.
#[derive(Debug)]
pub enum Polarization {
    Linear(Angle),
    Circular(Handedness),
    Elliptical(Complex<f64>, Complex<f64>),
}

/// The handedness of a circularly polarized beam.
///
/// A circularly polarized beam may either be left- or right-hand circularly polarized,
/// as determined by the right hand rule.
#[derive(Debug)]
pub enum Handedness {
    Left,
    Right,
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
#[derive(Debug, Clone)]
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
        let rel_phase = (y_phase - x_phase) % (2_f64 * pi);
        let y_rel = Complex::from_polar(&y_mag, &rel_phase);
        let x = Complex::<f64>::new(x_mag, 0.0);
        Beam {
            vec: ComplexVector::new(x, y_rel),
        }
    }

    fn remove_common_phase_mut(&mut self) {
        let (x_mag, x_phase) = self.vec.x.to_polar();
        let (y_mag, y_phase) = self.vec.y.to_polar();
        let rel_phase = (y_phase - x_phase) % (2_f64 * pi);
        self.vec.x = Complex::from_polar(&x_mag, &0_f64);
        self.vec.y = Complex::from_polar(&y_mag, &rel_phase);
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
prop_compose! {
    [pub(crate)] fn well_behaved_doubles()(x in (POSITIVE | NEGATIVE | ZERO).prop_filter(
            "Floats should be zero, or between 1e-12 and 1e12",
            |f| (f == &0_f64) || ((f.abs() < 1e12_f64) && (f.abs() > 1e-12_f64)))
    ) -> f64 {
        x
    }
}

#[cfg(test)]
prop_compose! {
    [pub(crate)] fn well_behaved_complexes()(x in well_behaved_doubles(),
                                y in well_behaved_doubles(),
    ) -> Complex<f64> {
        Complex::new(x, y)
    }
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

// Beam tests
#[cfg(test)]
proptest! {
   #[test]
   fn test_intensity_one_real_component(n in well_behaved_doubles()) {
       let beam = Beam::new(
           Complex::new(n, 0_f64),
           Complex::new(0_f64, 0_f64),
       );
       let intensity = beam.intensity().unwrap();
       assert_approx_eq!(n.powi(2), intensity);
   }

   #[test]
   fn test_intensity_is_never_negative(x in well_behaved_complexes(),
                                       y in well_behaved_complexes()) {
       let v = Beam::new(x, y);
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
       let beam = Beam::new(x, y);
       prop_assume!(beam.intensity().is_ok());
       assert_approx_eq!(beam.intensity().unwrap(), by_hand);
   }

   #[test]
   fn test_common_phase(x in well_behaved_complexes(),
                        y in well_behaved_complexes()) {
       let y_phase_new = y.arg() - x.arg();
       let beam = Beam::new(x, y);
       let new_beam = beam.remove_common_phase();
       assert_approx_eq!(new_beam.vec.x.arg(), 0_f64);
       assert_approx_eq!(new_beam.vec.y.arg(), y_phase_new % (2_f64 * pi));
   }

   #[test]
   fn test_common_phase_mut(x in well_behaved_complexes(),
                            y in well_behaved_complexes()) {
       let y_phase_new = y.arg() - x.arg();
       let mut beam = Beam::new(x, y);
       beam.remove_common_phase_mut();
       assert_approx_eq!(beam.vec.x.arg(), 0_f64);
       assert_approx_eq!(beam.vec.y.arg(), y_phase_new % (2_f64 * pi));
   }

   #[test]
   fn test_relative_phase(x in well_behaved_complexes(),
                          y in well_behaved_complexes()) {
       let expected_rad = y.arg() - x.arg();
       let expected_deg = expected_rad.to_degrees();
       let beam = Beam::new(x, y);
       assert_approx_eq!(expected_rad, beam.relative_phase_rad());
       assert_approx_eq!(expected_deg, beam.relative_phase_deg());
   }

    #[test]
    fn test_xy(x in well_behaved_complexes(),
               y in well_behaved_complexes()) {
        let beam = Beam::new(x, y);
        assert_eq!(beam.x(), x);
        assert_eq!(beam.y(), y);
    }
}

// JonesMatrix tests
#[cfg(test)]
proptest! {

   #[test]
   fn test_rotate_360_degrees_returns_original(m00 in well_behaved_complexes(),
                                               m01 in well_behaved_complexes(),
                                               m10 in well_behaved_complexes(),
                                               m11 in well_behaved_complexes()) {
       let mat = Matrix2::new(m00, m01, m10, m11);
       let angle = Angle::Degrees(360_f64);
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
       let angle = Angle::Degrees(0_f64);
       let rotated = rotate_matrix(&mat, &angle);
       assert_matrix_approx_eq!(mat, rotated);
   }

   #[test]
   fn test_rotate_0_rad_returns_original(m00 in well_behaved_complexes(),
                                         m01 in well_behaved_complexes(),
                                         m10 in well_behaved_complexes(),
                                         m11 in well_behaved_complexes()) {
       let mat = Matrix2::new(m00, m01, m10, m11);
       let angle = Angle::Radians(0_f64);
       let rotated = rotate_matrix(&mat, &angle);
       assert_matrix_approx_eq!(mat, rotated);
   }

}
