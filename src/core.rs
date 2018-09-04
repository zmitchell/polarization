#[derive(Debug)]
pub enum BeamPol {
    Linear(Angle),
    Circular(Handedness),
    Unpolarized,
}

#[derive(Debug)]
pub enum Handedness {
    Left,
    Right,
}

#[derive(Debug, Copy, Clone)]
pub enum Angle {
    Degrees(f64),
    Radians(f64),
}
