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
    use jones::common::{well_behaved_complexes, Beam, JonesVector};

    proptest! {
        #[test]
        fn test_identity_element_returns_beam(x in well_behaved_complexes(),
                                              y in well_behaved_complexes()) {
            let beam = Beam::new(x, y);
            let ident = IdentityElement::new();
            let beam_after = beam.apply_element(ident);
            assert_beam_approx_eq!(beam_after, beam);
        }

        #[test]
        fn test_identity_preserved_under_rotation(theta in 0_f64..90_f64) {
            let beam = Beam::linear(Angle::Degrees(theta));
            let ident = IdentityElement::new().rotated(Angle::Degrees(theta));
            let beam_after = beam.apply_element(ident);
            assert_beam_approx_eq!(beam_after, beam);
        }
    }
}
