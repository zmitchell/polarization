#[macro_use]
extern crate assert_approx_eq;
extern crate nalgebra as na;
#[macro_use]
extern crate derive_more;
extern crate num;
#[macro_use]
extern crate proptest;

mod core;
pub mod jones;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
