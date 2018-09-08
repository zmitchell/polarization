use na::{Matrix2, Vector2};
use num::complex::Complex;

use super::common::{
    rotate_matrix, well_behaved_complexes, well_behaved_doubles, Angle, ComplexMatrix,
    ElementParams, JonesError, JonesMatrix, JonesVector, MissingParameter, Result,
};

#[derive(Debug, Copy, Clone)]
pub(crate) struct CompositeElement {
    mat: ComplexMatrix,
}

impl CompositeElement {
    pub(crate) fn from_matrix(mat: ComplexMatrix) -> Self {
        CompositeElement { mat }
    }
}

impl JonesMatrix for CompositeElement {
    fn rotated(&self, angle: Angle) -> Self {
        CompositeElement {
            mat: rotate_matrix(&self.mat, &angle),
        }
    }

    fn rotate(&mut self, angle: Angle) {
        self.mat = rotate_matrix(&self.mat, &angle);
    }

    fn matrix(&self) -> ComplexMatrix {
        self.mat
    }
}
