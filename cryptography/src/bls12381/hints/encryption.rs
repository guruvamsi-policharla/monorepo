use aes_gcm::{aead::Aead, Aes256Gcm, Key, KeyInit};
use commonware_math::algebra::{CryptoGroup, Random};
use hkdf::Hkdf;
use rand::thread_rng;
use sha2::Sha256;

use crate::bls12381::{
    hints::{aggregate::EncryptionKey, crs::CRS, types::Ciphertext},
    primitives::{group::Scalar, variant::Variant},
};

/// t is the threshold for encryption and apk is the aggregated public key
pub fn encrypt<V: Variant>(
    ek: &EncryptionKey<V>,
    t: usize,
    crs: &CRS<V>,
    gamma_g2: V::Signature, // this should be hash_to_point(attestation_data)
    m: &[u8],
) -> Ciphertext<V> {
    let mut rng = &mut thread_rng();

    let g = crs.powers_of_g[0];
    let h = crs.powers_of_h[0];

    let mut sa1 = [V::Public::generator(); 2];
    let mut sa2 = [V::Signature::generator(); 6];

    let s = (0..5).map(|_| Scalar::random(&mut rng)).collect::<Vec<_>>();

    // sa1[0] = s0*ask + s3*g^{tau^{t}} + s4*g
    sa1[0] = (ek.ask * &s[0]) + &(crs.powers_of_g[t] * &s[3]) + &(crs.powers_of_g[0] * &s[4]);

    // sa1[1] = s2*g
    sa1[1] = g * &s[2];

    // sa2[0] = s0*h + s2*gamma_g2
    sa2[0] = (h * &s[0]) + &(gamma_g2 * &s[2]);

    // sa2[1] = s0*z_g2
    sa2[1] = ek.z_g2 * &s[0];

    // sa2[2] = s0*h^tau + s1*h^{tau^2}
    sa2[2] = crs.powers_of_h[1] * &s[0] + &(crs.powers_of_h[2] * &s[1]);

    // sa2[3] = s1*h
    sa2[3] = h * &s[1];

    // sa2[4] = s3*h
    sa2[4] = h * &s[3];

    // sa2[5] = s4*h^{tau}
    sa2[5] = (crs.powers_of_h[1]) * &s[4];

    // NOTE: BLST doesn't implement GT operations
    // So we do it suboptimally via a pairing
    // enc_key = s4*e_gh
    // let enc_key = ek.e_gh * s[4];
    let enc_key = V::pairing(
        &(V::Public::generator() * &s[4]),
        &V::Signature::generator(),
    );
    let enc_key_bytes = enc_key.as_slice();

    // derive an encapsulation key from enc_key using an HKDF
    let hk = Hkdf::<Sha256>::new(None, &enc_key_bytes);
    let mut aes_key = [0u8; 32];
    let mut aes_nonce = [0u8; 12];
    hk.expand(&[1], &mut aes_key).unwrap();
    hk.expand(&[2], &mut aes_nonce).unwrap();

    // encrypt the message m using the derived key
    let aes_key: &Key<Aes256Gcm> = &aes_key.into();
    let cipher = Aes256Gcm::new(aes_key);
    let ct = cipher.encrypt(&aes_nonce.into(), m).unwrap();

    Ciphertext {
        gamma_g2,
        sa1,
        sa2,
        ct,
        t,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bls12381::{
        hints::{
            aggregate::AggregateKey,
            setup::{LagPublicKey, SecretKey},
        },
        primitives::variant::{self, MinPk},
    };

    #[test]
    fn test_encryption() {
        let mut rng = thread_rng();
        let n = 8;
        let crs = CRS::new(n);

        let msg = b"Hello, world!";

        let mut sk: Vec<SecretKey> = Vec::new();
        let mut pk: Vec<LagPublicKey<MinPk>> = Vec::new();

        for i in 0..n {
            sk.push(SecretKey::new(i));
            pk.push(sk[i].get_lagrange_pk(i, &crs))
        }

        let (_ak, ek) = AggregateKey::new(pk, &crs);

        // todo: should be hash to curve
        let gamma_g2 =
            <MinPk as variant::Variant>::Signature::generator() * &Scalar::random(&mut rng);

        let _ct = encrypt::<MinPk>(&ek, 2, &crs, gamma_g2, msg);

        // let mut ct_bytes = Vec::new();
        // ct.serialize_compressed(&mut ct_bytes).unwrap();
        // println!("Compressed ciphertext: {} bytes", ct_bytes.len());

        // let mut g1_bytes = Vec::new();
        // let mut g2_bytes = Vec::new();
        // let mut e_gh_bytes = Vec::new();

        // let g = G1::generator();
        // let h = G2::generator();

        // g.serialize_compressed(&mut g1_bytes).unwrap();
        // h.serialize_compressed(&mut g2_bytes).unwrap();
        // ek.e_gh.serialize_compressed(&mut e_gh_bytes).unwrap();

        // println!("G1 len: {} bytes", g1_bytes.len());
        // println!("G2 len: {} bytes", g2_bytes.len());
        // println!("GT len: {} bytes", e_gh_bytes.len());
    }
}
