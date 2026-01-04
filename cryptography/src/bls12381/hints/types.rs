//! common types
use crate::bls12381::primitives::variant::Variant;

#[derive(Clone, PartialEq)]
pub struct Ciphertext<V: Variant> {
    pub gamma_g2: V::Signature,
    pub sa1: [V::Public; 2],
    pub sa2: [V::Signature; 6],
    pub ct: Vec<u8>, //encrypted message
    pub t: usize,    //threshold
}

impl<V: Variant> Ciphertext<V> {
    pub fn new(
        gamma_g2: V::Signature,
        sa1: [V::Public; 2],
        sa2: [V::Signature; 6],
        ct: Vec<u8>,
        t: usize,
    ) -> Self {
        Ciphertext {
            gamma_g2,
            sa1,
            sa2,
            ct,
            t,
        }
    }
}
