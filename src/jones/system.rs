//! An optical system that encapsulates a beam and the elements that it will pass
//! through.

use super::common::{Beam, JonesError, JonesMatrix, JonesVector, Result};
use super::composite::CompositeElement;
use super::hwp::HalfWavePlate;
use super::identity::IdentityElement;
use super::polarizer::Polarizer;
use super::qwp::QuarterWavePlate;
use super::retarder::Retarder;
use super::rotator::PolarizationRotator;
#[cfg(test)]
use proptest::option;
#[cfg(test)]
use proptest::prelude::*;

/// A type that represents the various optical elements that may appear in the
/// optical system.
#[derive(Debug, Copy, Clone)]
pub enum OpticalElement {
    /// An ideal linear polarizer.
    Polarizer(Polarizer),
    /// An optical element that rotates the polarization of a beam.
    PolarizationRotator(PolarizationRotator),
    /// An optical retarder with an arbitrary phase delay.
    Retarder(Retarder),
    /// An ideal quarter-wave plate.
    QuarterWavePlate(QuarterWavePlate),
    /// An ideal half-wave plate.
    HalfWavePlate(HalfWavePlate),
    /// An element that passes a beam through untouched.
    Identity(IdentityElement),
    /// An element that represents the composition of several other elements.
    Composite(CompositeElement),
}

#[cfg(test)]
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum OpticalElementType {
    Polarizer,
    PolarizationRotator,
    Retarder,
    QuarterWavePlate,
    HalfWavePlate,
    Identity,
    Composite,
    Any,
}

#[cfg(test)]
impl Default for OpticalElementType {
    fn default() -> Self {
        OpticalElementType::Any
    }
}

#[cfg(test)]
impl Arbitrary for OpticalElement {
    type Parameters = OpticalElementType;
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(args: Self::Parameters) -> Self::Strategy {
        match args {
            OpticalElementType::Polarizer => any::<Polarizer>()
                .prop_map(|x| OpticalElement::Polarizer(x))
                .boxed(),
            OpticalElementType::PolarizationRotator => any::<PolarizationRotator>()
                .prop_map(|x| OpticalElement::PolarizationRotator(x))
                .boxed(),
            OpticalElementType::Retarder => any::<Retarder>()
                .prop_map(|x| OpticalElement::Retarder(x))
                .boxed(),
            OpticalElementType::QuarterWavePlate => any::<QuarterWavePlate>()
                .prop_map(|x| OpticalElement::QuarterWavePlate(x))
                .boxed(),
            OpticalElementType::HalfWavePlate => any::<HalfWavePlate>()
                .prop_map(|x| OpticalElement::HalfWavePlate(x))
                .boxed(),
            OpticalElementType::Identity => any::<IdentityElement>()
                .prop_map(|x| OpticalElement::Identity(x))
                .boxed(),
            OpticalElementType::Composite => any::<CompositeElement>()
                .prop_map(|x| OpticalElement::Composite(x))
                .boxed(),
            OpticalElementType::Any => prop_oneof![
                any::<Polarizer>().prop_map(|x| OpticalElement::Polarizer(x)),
                any::<PolarizationRotator>().prop_map(|x| OpticalElement::PolarizationRotator(x)),
                any::<Retarder>().prop_map(|x| OpticalElement::Retarder(x)),
                any::<QuarterWavePlate>().prop_map(|x| OpticalElement::QuarterWavePlate(x)),
                any::<HalfWavePlate>().prop_map(|x| OpticalElement::HalfWavePlate(x)),
                any::<IdentityElement>().prop_map(|x| OpticalElement::Identity(x)),
                any::<CompositeElement>().prop_map(|x| OpticalElement::Composite(x)),
            ]
            .boxed(),
        }
        .boxed()
    }
}

/// A type that contains a beam and the elements that it will pass through.
///
/// # Examples
/// An optical system is constructed using the builder pattern.
/// ```
/// # extern crate polarization;
/// use polarization::jones::*;
/// // A linearly polarized beam parallel to the x-axis
/// let beam = Beam::linear(Angle::Degrees(0.0));
/// // A linear polarizer at 45 degrees counter-clockwise from the x-axis
/// let pol = OpticalElement::Polarizer(Polarizer::new(Angle::Degrees(45.0)));
/// // The optical system
/// let system = OpticalSystem::new()
///     .with_beam(beam)
///     .with_element(pol);
/// ```
/// If you have several elements, you may also add them to the system using the
/// `with_elements` method.
///
/// The beam may be propagated through the elements in the system, returning a new beam.
/// ```
/// # use polarization::jones::*;
/// # let beam = Beam::linear(Angle::Degrees(0.0));
/// # let pol = OpticalElement::Polarizer(Polarizer::new(Angle::Degrees(45.0)));
/// # let system = OpticalSystem::new()
/// #     .with_beam(beam)
/// #     .with_element(pol);
/// let final_beam: Result<Beam> = system.propagate();
/// ```
/// This operation may fail because you're human and maybe you forgot to add a beam or
/// elements to your system before calling `propagate`.
///
/// The order in which you add elements to the system is the same order that the beam
/// will travel through the elements. For example, consider a horizontal beam, a
/// vertical polarizer, and a rotator that will rotate the beam by 90 degrees. If the
/// beam passes through the polarizer first, the intensity after the polarizer will be
/// zero since the polarization of the beam is perpendicular to the axis of the
/// polarizer.
/// ```
/// # use polarization::jones::*;
/// let beam = Beam::linear(Angle::Degrees(0.0));
/// let pol = OpticalElement::Polarizer(
///     Polarizer::new(Angle::Degrees(90.0))
/// );
/// let rot = OpticalElement::PolarizationRotator(
///     PolarizationRotator::new(Angle::Degrees(90.0))
/// );
/// let system = OpticalSystem::new()
///     .with_beam(beam.clone())
///     .with_element(pol.clone())  // this kills the beam!
///     .with_element(rot.clone());  // this does nothing now
/// let beam_after = system.propagate().unwrap();
/// assert!(beam_after.intensity().unwrap() < 1e-6);
/// ```
/// However, if the beam passes through the rotator first, the beam after the
/// rotator will be parallel to the polarizer, and will pass through unaffected.
/// ```
/// # use polarization::jones::*;
/// # let beam = Beam::linear(Angle::Degrees(0.0));
/// # let pol = OpticalElement::Polarizer(
/// #     Polarizer::new(Angle::Degrees(90.0))
/// # );
/// # let rot = OpticalElement::PolarizationRotator(
/// #     PolarizationRotator::new(Angle::Degrees(90.0))
/// # );
/// let system = OpticalSystem::new()
///     .with_beam(beam.clone())
///     .with_element(rot.clone())  // beam and polarizer are now parallel
///     .with_element(pol.clone());
/// let beam_after = system.propagate().unwrap();
/// assert!(beam_after.intensity().unwrap() > 0.99);
/// ```
#[derive(Debug, Clone)]
pub struct OpticalSystem {
    pub beam: Option<Beam>,
    pub elements: Option<Vec<OpticalElement>>,
}

impl OpticalSystem {
    /// Constructs a new optical system.
    ///
    /// The optical system constructed by this method does not contain a beam, and does not
    /// contain any optical elements. A beam and elements are both necessary for doing a
    /// polarization simulation.
    pub fn new() -> Self {
        OpticalSystem {
            beam: None,
            elements: None,
        }
    }

    /// Add a beam to the optical system.
    ///
    /// Calling this method more than once will overwrite the previously added beam.
    pub fn with_beam(self, beam: Beam) -> Self {
        OpticalSystem {
            beam: Some(beam),
            elements: self.elements,
        }
    }

    /// Add an optical element to the system.
    ///
    /// The beam will pass through the elements in the order that they are added.
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

    /// Add several optical elements to the system at once.
    ///
    /// If elements have already been added to the system, these elements will be added after
    /// the existing elements.
    pub fn with_elements(self, elems: Vec<OpticalElement>) -> Self {
        if self.elements.is_none() {
            let system = OpticalSystem {
                beam: self.beam,
                elements: Some(elems),
            };
            return system;
        }
        let mut elements = self.elements.clone().unwrap();
        elements.extend(elems);
        OpticalSystem {
            beam: self.beam,
            elements: Some(elements),
        }
    }

    /// Propagate the beam through the elements in the system, returning the final beam.
    /// Will return an error if there is no beam in the system, or if there are no elements in
    /// the system.
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
impl Arbitrary for OpticalSystem {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_: Self::Parameters) -> Self::Strategy {
        (
            option::of(any::<Beam>()),
            option::of(any::<Vec<OpticalElement>>()),
        )
            .prop_map(|(beam, elements)| OpticalSystem { beam, elements })
            .boxed()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use jones::common::{any_complex, Angle};
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
        fn test_beam_passes_through(x in any_complex(), y in any_complex(),) {
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
