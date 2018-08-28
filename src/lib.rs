#[macro_use] extern crate assert_approx_eq;
extern crate nalgebra as na;
#[macro_use] extern crate derive_more;
extern crate num;
#[macro_use] extern crate proptest;

pub mod jones;
mod core;


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
