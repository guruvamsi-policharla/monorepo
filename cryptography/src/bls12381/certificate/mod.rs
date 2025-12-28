//! BLS12-381 certificate scheme implementations.
//!
//! This module provides two signing scheme implementations:
//!
//! - [`multisig`]: Multi-signature scheme with attributable signatures
//! - [`threshold`]: Threshold signature scheme with non-attributable signatures
//! - [`silent`]: Threshold signature scheme with silent setup and optionally attributable signatures

pub mod multisig;
pub mod silent;
pub mod threshold;
