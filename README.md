# Polarization
[![Documentation](https://docs.rs/polarization/badge.svg)](https://docs.rs/polarization)
[![Crates.io](https://img.shields.io/crates/v/polarization.svg)](https://crates.io/polarization)
![Licenses](https://img.shields.io/crates/l/polarization.svg)

Have you ever wondered what would happen if you passed a linearly polarized beam
through a quarter-wave plate at 46 degrees rather than 45 degrees relative to the
fast axis of a quarter-wave plate? Who am I kidding, of course you have! This
library lets you pass a beam through several optical elements and see what comes
out the other side.

The canonical methods for simulating the polarization of a beam are
[Jones calculus](https://en.wikipedia.org/wiki/Jones_calculus) and
[Mueller calculus](https://en.wikipedia.org/wiki/Mueller_calculus), but only Jones calculus
is implemented at this point.

Currently there are several standard optical elements implemented, with support for reflections
from surfaces (dielectric and metallic) coming in the near future.
* Linear polarizer
* Polarization rotator
* Quarter-wave plate
* Half-wave plate
* Retarder

There is support for linearly polarized, circularly polarized, and arbitrarily
(elliptically) polarized beams.

For more details, check out the [documentation](https://docs.rs/polarization).

## Example
```rust
let beam = Beam::linear(Angle::Degrees(0.0));
let pol = OpticalElement::Polarizer(Polarizer::new(Angle::Degrees(45.0)));
let system = OpticalSystem::new()
    .with_beam(beam)
    .with_element(pol);
let final_beam: Result<Beam> = system.propagate();
let final_intensity: Result<f64> = final_beam.intensity();
```

## License

Licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.
