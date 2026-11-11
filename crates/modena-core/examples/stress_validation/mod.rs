//! Independent bounded support for tensor-pressure validation.
// Each integration test links only its own slice of this shared example module.
#![allow(dead_code)]
pub mod ase;
pub mod costs;
pub mod ensemble;
pub mod fixtures;
pub mod matched;
pub mod mesh_oracle;
pub mod ou;
pub mod protocol;
pub mod statistics;
