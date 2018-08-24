use na::{Matrix2, Vector2};
use num::complex::Complex;

use core::{Angle, BeamPol, Handedness};


pub trait JonesVector {
    /// The intensity of the beam represented by the Jones vector. Note that
    /// this is simply `V* x V`.
    fn intensity(&self) -> f64;

    /// For the vector `V = (A, B)`, return a new vector with the phase common
    /// to `A` and `B` removed. For the returned vector `V' = (A', B')`, `A'`
    /// will be real, and `B'` will have some phase relative to `A'`.
    fn remove_common_phase(&self) -> Self;

    /// For the vector `V = (A, B)`, remove the phase common to both `A` and `B`
    /// in place. See `remove_common_phase`.
    fn remove_common_phase_mut(&mut self);

    /// Returns the relative phase between the x- and y-components of the vector
    /// in degrees.
    fn relative_phase_deg(&self) -> Angle;

    /// Returns the relative phase between the x- and y-components of the vector
    /// in radians.
    fn relative_phase_rad(&self) -> Angle;
}

impl JonesVector for Vector2<Complex<f64>> {
    fn intensity(&self) -> f64 {
        let conj = Vector2::new(
            self.x.conj(),
            self.y.conj(),
        );
        conj.dot(self).norm_sqr().sqrt()
    }

    fn remove_common_phase(&self) -> Self {
        let (x_mag, x_phase) = self.x.to_polar();
        let (y_mag, y_phase) = self.y.to_polar();
        let rel_phase = y_phase - x_phase;
        let y_rel = Complex::from_polar(&y_mag, &rel_phase);
        Vector2::new(
            Complex::new(x_mag, 0.0),
            y_rel,
        )
    }

    fn remove_common_phase_mut(&mut self) {
        let (x_mag, x_phase) = self.x.to_polar();
        let (y_mag, y_phase) = self.y.to_polar();
        let rel_phase = y_phase - x_phase;
        self.x = Complex::from_polar(&x_mag, &0.0);
        self.y = Complex::from_polar(&y_mag, &rel_phase);
    }

    fn relative_phase_deg(&self) -> Angle {
        let (_, x_phase) = self.x.to_polar();
        let (_, y_phase) = self.y.to_polar();
        let phase = y_phase.to_degrees() - x_phase.to_degrees();
        Angle::Degrees(phase)
    }

    fn relative_phase_rad(&self) -> Angle {
        let (_, x_phase) = self.x.to_polar();
        let (_, y_phase) = self.y.to_polar();
        let phase = y_phase - x_phase;
        Angle::Radians(phase)
    }
}

/// Produces the Jones vector for a linearly polarized beam that
/// is parallel to the x-axis.
pub fn beam_horizontal() -> Vector2<Complex<f64>> {
    Vector2::new(
        Complex::new(1.0, 0.0),
        Complex::new(0.0, 0.0),
    )
}

/// Produces the Jones vector for a linearly polarized beam that
/// is perpendicular to the x-axis.
pub fn beam_vertical() -> Vector2<Complex<f64>> {
    Vector2::new(
        Complex::new(0.0, 0.0),
        Complex::new(1.0, 0.0),
    )
}

