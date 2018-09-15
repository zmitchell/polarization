//! Polarization simulations using [Jones calculus](https://en.wikipedia.org/wiki/Jones_calculus).
//!
//! Jones calculus represents the components of an electric field (its polarization) as
//! a 2x1 column vector (a Jones vector), and the effects of optical elements as 2x2
//! matrices (Jones matrices). The final polarization is found by multiplying the Jones
//! vector of the beam by the Jones matrices of all the optical elements that the beam
//! passes through. In general, both Jones vectors and Jones matrices are complex (have
//! imaginary components), where the imaginary parts encode the phase of each component.
//!
//! # Conventions
//! This is physics, so of course there are going to be some conventions.
//!
//! ### Axes
//! The beam propagates in the +z-direction. This means that the polarization of the
//! beam lives in the xy-plane.
//!
//! ### Angles
//! All angles are measured counter-clockwise from the x-axis.
//!
//! # Overview
//! A basic polarization simulation involves just a few steps:
//! * Create a beam
//! * Create some optical elements
//! * Create an `OpticalSystem`
//! * Add the beam and the elements to the `OpticalSystem`
//! * Call `propagate` to send the beam through the elements
//!
//! ### Elements
//! There are several standard optical elements defined by the `jones` module:
//! * linear polarizer (`Polarizer`)
//! * half-wave plate (`HalfWavePlate`)
//! * quarter-wave plate (`QuarterWavePlate`)
//! * retarder (`Retarder`)
//! * polarization rotator (`PolarizationRotator`)
//!
//! Support for reflections from a dielectric or metal surface will be added in the near
//! future.
//!
//! What if you want to use an optical element type that isn't included in this library?
//! There is a slight hack that will let you do this. There is an element type called
//! `CompositeElement` that is used to represent the effects of several elements with a
//! single Jones matrix. This element type is constructed from a user-supplied matrix.
//! If you can create the matrix for your custom element, then you can wrap it in a
//! `CompositeElement` and use it in a simulation.
//!
//! ### Beams
//! A beam can be constructed in a number of ways. There are convenience methods for
//! defining common polarization states:
//! * linearly polarized -> `Beam::linear`
//! * circularly polarized -> `Beam::circular`
//!
//! You may also construct an arbitrary beam using `Beam::new`.
//!
//! ### Optical System
//! An optical system is the construct that performs the simulation. To perform a
//! simulation, you must add a beam and some elements to the optical system, then call
//! `propagate`.
//! ```
//! # use polarization::jones::*;
//! let beam = Beam::linear(Angle::Degrees(0.0));
//! let pol = OpticalElement::Polarizer(Polarizer::new(Angle::Degrees(45.0)));
//! let system = OpticalSystem::new()
//!     .with_beam(beam)
//!     .with_element(pol);
//! let final_beam: Result<Beam> = system.propagate();
//! ```
//!
//! If you add multiple elements to the optical system, the beam will travel through
//! them in the order that they are added to the optical system.

pub mod common;
pub mod composite;
pub mod hwp;
pub mod identity;
pub mod polarizer;
pub mod qwp;
pub mod retarder;
pub mod rotator;
pub mod system;

pub use self::common::*;
pub use self::composite::*;
pub use self::hwp::*;
pub use self::identity::*;
pub use self::polarizer::*;
pub use self::qwp::*;
pub use self::retarder::*;
pub use self::rotator::*;
pub use self::system::*;
