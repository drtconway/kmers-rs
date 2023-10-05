#![warn(missing_docs)]

//! K-mers and associated operations.
//! 
//! This library provides functionality for extracting k-mers from sequences,
//! and manipulating them in useful ways. The underlying representation is
//! 64-bit integers (`u64`), so k > 32 is not supported by this library.
//! 
//! K-mers (or q-grams in some computer science contexts) are k-length sequences
//! of DNA/RNA "letters" represented as unsigned integers. Following usual practice,
//! 
//! * "A" -> b00
//! * "C" -> b01
//! * "G" -> b10
//! * "T" or "U" -> b11
//! 
//! which has the nice property that the complementary bases are bitwise complements.
//!  

mod basics;
mod accumulator;
mod frequency;
mod dot;
mod jaccard;

pub use basics::*;
pub use accumulator::*;
pub use frequency::*;
pub use dot::*;
pub use jaccard::*;