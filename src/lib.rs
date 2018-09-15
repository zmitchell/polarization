//! Simulate the polarization of a laser beam.
//!
//! Have you ever wondered what would happen if you passed a linearly polarized beam
//! through a quarter-wave plate at 46 degrees rather than 45 degrees relative to the
//! fast axis of a quarter-wave plate? Who am I kidding, of course you have! This
//! library lets you pass a beam through several optical elements and see what comes
//! out the other side.
//!
//! Currently there is only support for
//! [Jones calculus](https://en.wikipedia.org/wiki/Jones_calculus), but support for
//! [Mueller calculus](https://en.wikipedia.org/wiki/Mueller_calculus) may be added at
//! some point in the future.
//!
//! For an overview of what you can do with this crate, check out the `jones` module.
#![cfg_attr(feature = "cargo-clippy", allow(needless_pass_by_value))]
#![allow(dead_code)]
#[cfg(test)]
#[macro_use]
extern crate assert_approx_eq;
extern crate nalgebra as na;
extern crate num;
#[cfg(test)]
#[macro_use]
extern crate proptest;

pub mod jones;
