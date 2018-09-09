#![cfg_attr(feature = "cargo-clippy", allow(needless_pass_by_value))]
#[cfg(test)]
#[macro_use]
extern crate assert_approx_eq;
extern crate nalgebra as na;
extern crate num;
#[cfg(test)]
#[macro_use]
extern crate proptest;

pub mod jones;
