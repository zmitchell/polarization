use super::common::{Beam, JonesError, JonesMatrix, JonesVector, Result};
use super::composite::CompositeElement;
use super::hwp::HalfWavePlate;
use super::identity::IdentityElement;
use super::polarizer::Polarizer;
use super::qwp::QuarterWavePlate;
use super::retarder::Retarder;
use super::rotator::PolarizationRotator;

#[derive(Debug, Copy, Clone)]
pub enum OpticalElement {
    Polarizer(Polarizer),
    PolarizationRotator(PolarizationRotator),
    Retarder(Retarder),
    QuarterWavePlate(QuarterWavePlate),
    HalfWavePlate(HalfWavePlate),
    Identity(IdentityElement),
    Composite(CompositeElement),
}

#[derive(Debug, Clone)]
pub struct OpticalSystem {
    pub beam: Option<Beam>,
    pub elements: Option<Vec<OpticalElement>>,
}

impl OpticalSystem {
    pub fn new() -> Self {
        OpticalSystem {
            beam: None,
            elements: None,
        }
    }

    pub fn with_beam(self, beam: Beam) -> Self {
        OpticalSystem {
            beam: Some(beam),
            elements: self.elements,
        }
    }

    pub fn with_element(self, elem: OpticalElement) -> Self {
        let elements = match self.elements {
            Some(elements) => {
                elements.clone().push(elem);
                elements
            }
            None => vec![elem],
        };
        OpticalSystem {
            beam: self.beam,
            elements: Some(elements),
        }
    }

    /// Propagate the beam through the elements in the system, returning the final beam.
    /// Will return an error if there is no beam in the system, or if there are no elements in the system.
    pub fn propagate(&self) -> Result<Beam> {
        // Bring the variant names into scope just for convenience.
        use self::OpticalElement::*;
        if self.elements.is_none() {
            return Err(JonesError::NoElements);
        }
        if self.beam.is_none() {
            return Err(JonesError::NoBeam);
        }
        // This will be the "accumulator" for the call to 'fold'
        let ident = IdentityElement::new();
        // You can combine the optical elements into a single element by multiplying
        // them all together. Then the new beam can be created with 'apply_element'
        let composite_mat =
            self.elements
                .clone()
                .unwrap()
                .iter()
                .fold(ident.matrix(), |acc, elem| {
                    let mat = match *elem {
                        Polarizer(pol) => pol.matrix(),
                        PolarizationRotator(rot) => rot.matrix(),
                        Retarder(ret) => ret.matrix(),
                        QuarterWavePlate(qwp) => qwp.matrix(),
                        HalfWavePlate(hwp) => hwp.matrix(),
                        Identity(id) => id.matrix(),
                        Composite(comp) => comp.matrix(),
                    };
                    mat * acc
                });
        let composite = CompositeElement::from_matrix(composite_mat);
        Ok(self.beam.clone().unwrap().apply_element(composite))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[macro_use]
    use jones::common::{Angle, well_behaved_complexes};
    use num::complex::Complex;

    #[test]
    fn test_elements_applied_in_order() {
        // Pass a horizontal beam through a vertical polarizer and then a polarization
        // rotator. If applied in the correct order, the polarizer will kill the beam
        // since the beam and the polarizer axis are perpendicular. The polarization
        // rotator will have no effect in this case. If the elements are applied in the
        // wrong order, then the polarization rotator will rotate the beam to be
        // parallel to the polarizer and the final beam will have a non-zero intensity.
        let beam = Beam::linear(Angle::Degrees(0.0));
        let zero_beam = Beam::new(Complex::new(0_f64, 0_f64), Complex::new(0_f64, 0_f64));
        let pol = Polarizer::new(Angle::Degrees(90.0));
        let rot = PolarizationRotator::new(Angle::Degrees(90.0));
        let system = OpticalSystem::new()
            .with_beam(beam)
            .with_element(OpticalElement::Polarizer(pol))
            .with_element(OpticalElement::PolarizationRotator(rot));
        let beam_after = system.propagate();
        assert!(beam_after.is_ok());
        assert_beam_approx_eq!(beam_after.unwrap(), zero_beam);
    }

    proptest!{
        #[test]
        fn test_beam_passes_through(x in well_behaved_complexes(),
                                    y in well_behaved_complexes(),) {
            let beam = Beam::new(x, y);
            let ident = IdentityElement::new();
            let element = OpticalElement::Identity(ident);
            let system = OpticalSystem::new()
                .with_beam(beam.clone())
                .with_element(element);
            let beam_after = system.propagate();
            assert!(beam_after.is_ok());
            assert_beam_approx_eq!(beam_after.unwrap(), beam);
        }

    }
}
