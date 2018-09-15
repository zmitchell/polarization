//! Implementation of an element that composes several other elements.
//!
//! One of the conveniences of Jones calculus is the ability to compose several optical
//! elements into a single element by multiplying their matrices together. This element
//! represents the result of composing several elements together, and, as a result, may
//! only be constructed from an existing `ComplexMatrix`.

use super::common::{rotate_matrix, Angle, ComplexMatrix, JonesMatrix};

/// An optical element that represents the composition of several other elements.
///
/// See the module-level documentation for more information.
#[derive(Debug, Copy, Clone)]
pub struct CompositeElement {
    mat: ComplexMatrix,
}

impl CompositeElement {
    /// Constructs a new `CompositeElement` from an existing `ComplexMatrix`.
    ///
    /// It is assumed that this matrix is the result of multiplying other Jones matrices
    /// together. No checks are performed on `mat` though, so you may supply an
    /// arbitrary Jones matrix to define your own optical element.
    pub fn from_matrix(mat: ComplexMatrix) -> Self {
        CompositeElement { mat }
    }
}

impl JonesMatrix for CompositeElement {
    /// Returns the element rotated counter-clockwise by `angle`.
    fn rotated(&self, angle: Angle) -> Self {
        CompositeElement {
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
