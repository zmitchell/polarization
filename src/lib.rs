extern crate nalgebra as na;
#[macro_use] extern crate derive_more;
extern crate num;

pub mod jones;
mod core;


#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
