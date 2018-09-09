use super::common::{rotate_matrix, Angle, ComplexMatrix, JonesMatrix};

#[derive(Debug, Copy, Clone)]
pub struct CompositeElement {
    mat: ComplexMatrix,
}

impl CompositeElement {
    pub fn from_matrix(mat: ComplexMatrix) -> Self {
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
