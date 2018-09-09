use num::complex::Complex;

use super::common::{rotate_matrix, Angle, ComplexMatrix, JonesMatrix};

#[derive(Debug, Copy, Clone)]
pub struct IdentityElement {
    mat: ComplexMatrix,
}

impl IdentityElement {
    pub fn new() -> Self {
        let zero = Complex::new(0_f64, 0_f64);
        let one = Complex::new(1_f64, 0_f64);
        let mat = ComplexMatrix::new(one, zero, zero, one);
        IdentityElement { mat }
    }
}

impl JonesMatrix for IdentityElement {
    fn rotated(&self, angle: Angle) -> Self {
        let mat = rotate_matrix(&self.matrix(), &angle);
        IdentityElement { mat }
    }

    fn rotate(&mut self, angle: Angle) {
        self.mat = rotate_matrix(&self.matrix(), &angle);
    }

    fn matrix(&self) -> ComplexMatrix {
        self.mat
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[macro_use]
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
