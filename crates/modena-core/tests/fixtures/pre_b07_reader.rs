//! Compile against the unmodified d6cd33ca baseline modena-core.
//! Arguments are a newly written v2 rotor checkpoint and a scratch output path.
use modena_core::io::restart::{load, save};

fn main() {
    let paths: Vec<_> = std::env::args_os().skip(1).collect();
    let cp = load(&paths[0]).expect("pre-B07 reader must accept the new v2 writer");
    let state = cp.continuation.as_ref().unwrap();
    assert_eq!((cp.step, cp.atoms, state.constraint_count), (17, 2, 1));
    assert_eq!(state.integrator["kind"], "langevin_baoab");
    save(&paths[1], &cp).unwrap();
    let bytes = std::fs::read(&paths[0]).unwrap();
    assert_eq!(bytes, std::fs::read(&paths[1]).unwrap());
    println!("original v2 reader/writer round-tripped {} bytes exactly, step {}, atoms {}, RNG word position {}", bytes.len(), cp.step, cp.atoms, state.integrator["rng"]["word_pos"]);
}
