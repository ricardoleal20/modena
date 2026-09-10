//! Error types for the MODENA core.

use thiserror::Error;

#[derive(Error, Debug)]
pub enum ModenaError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("TOML parse error: {0}")]
    Toml(#[from] toml::de::Error),

    #[error("configuration error: {0}")]
    Config(String),

    #[error("invalid input file: {0}")]
    Input(String),

    #[error("physics error: {0}")]
    Physics(String),

    #[error("atom count mismatch: expected {expected}, found {found}")]
    AtomCount { expected: usize, found: usize },
}

pub type Result<T> = std::result::Result<T, ModenaError>;
