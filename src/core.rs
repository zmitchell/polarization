pub enum BeamPol {
    Linear(Angle),
    Circular(Handedness),
    Unpolarized,
}

pub enum Handedness {
    Left,
    Right,
}

pub enum Angle {
    Degrees(f64),
    Radians(f64),
}
