//! An element that has no effect.
//!
//! An identity element simply passes the incident beam straight through. In real life
//! this could simply be the beam traveling through a transparent medium. This element
//! is mostly useful for testing and other internal uses in `polarization`, but may
//! prove handy for other uses as well.

use num::complex::Complex;

use super::common::{rotate_matrix, Angle, ComplexMatrix, JonesMatrix};
#[cfg(test)]
use proptest::prelude::*;

/// An optical element that simply passes a beam through untouched.
#[derive(Debug, Copy, Clone)]
pub struct IdentityElement {
    mat: ComplexMatrix,
}

impl IdentityElement {
    /// Constructs a new `IdentityElement`.
    ///
    /// There are no parameters because this element doesn't really do anything.
    pub fn new() -> Self {
        let zero = Complex::new(0_f64, 0_f64);
        let one = Complex::new(1_f64, 0_f64);
        let mat = ComplexMatrix::new(one, zero, zero, one);
        IdentityElement { mat }
    }
}

impl JonesMatrix for IdentityElement {
    /// Returns an `IdentityElement` that has been rotated counter-clockwise by
    /// `angle`.
    ///
    /// Note that the identity element only passes the original beam through when the
    /// element and the beam are in the same coordinate system. For example, if you
    /// start with a given beam and pass it through an identity element that has
    /// been rotated, in general you will not get back the same beam. However, if you
    /// start with a given beam, rotate both the beam and the element by the same angle,
    /// then rotate the beam back to its original orientation, you will recover the beam
    /// that you started with.
    fn rotated(&self, angle: Angle) -> Self {
        let mat = rotate_matrix(&self.matrix(), &angle);
        IdentityElement { mat }
    }

    /// Replaces this element with one that has been rotated counter-clockwise by
    /// `angle`.
    ///
    /// See above for notes on the semantics of rotating an `IdentityElement`.
    fn rotate(&mut self, angle: Angle) {
        self.mat = rotate_matrix(&self.matrix(), &angle);
    }

    /// Returns the 2x2 Jones matrix of the element.
    fn matrix(&self) -> ComplexMatrix {
        self.mat
    }
}

#[cfg(test)]
impl Arbitrary for IdentityElement {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        Just(IdentityElement::new()).boxed()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use jones::common::{float_angle, Beam, JonesVector};

    proptest! {
        #[test]
        fn test_identity_element_returns_same_beam(beam: Beam) {
            let ident = IdentityElement::new();
            let beam_after = beam.apply_element(ident);
            prop_assert_beam_approx_eq!(beam_after, beam);
        }

        #[test]
        fn test_multiple_identity_elements_preserves_beam(beam: Beam, elems: Vec<IdentityElement>) {
            let mut beam_after = beam.clone();
            for elem in elems {
                beam_after = beam_after.apply_element(elem);
            }
            prop_assert_beam_approx_eq!(beam_after, beam);
        }

        #[test]
        fn test_identity_preserved_under_rotation(theta in float_angle()) {
            let beam = Beam::linear(Angle::Degrees(theta));
            let ident = IdentityElement::new().rotated(Angle::Degrees(theta));
            let beam_after = beam.apply_element(ident);
            prop_assert_beam_approx_eq!(beam_after, beam);
        }

    }
}
