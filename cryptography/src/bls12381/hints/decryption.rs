use aes_gcm::{aead::Aead, Aes256Gcm, Key, KeyInit};
use commonware_math::{
    algebra::{Additive, Field, Ring, Space},
    poly::Poly,
};
use commonware_utils::vec::NonEmptyVec;
use hkdf::Hkdf;
use sha2::Sha256;

use crate::bls12381::{
    hints::{
        aggregate::AggregateKey, crs::CRS, fft_settings::Settings, setup::PartialDecryption,
        types::Ciphertext, utils::interp_mostly_zero,
    },
    primitives::{group::Scalar, variant::Variant},
};

pub fn agg_dec<V: Variant>(
    partial_decryptions: &[PartialDecryption<V>], /* insert 0 if a party did not respond or
                                                   * verification failed */
    ct: &Ciphertext<V>,
    selector: &[bool],
    agg_key: &AggregateKey<V>,
    crs: &CRS<V>,
) -> Vec<u8> {
    let domain = Settings::new(crs.n.trailing_zeros() as usize).unwrap();
    // Radix2EvaluationDomain::<E::ScalarField>::new(crs.n).unwrap();
    // let domain_elements: Vec<E::ScalarField> = domain.elements().collect();

    // points is where B is set to zero
    // parties is the set of parties who have signed
    let mut points = vec![];
    let mut parties: Vec<usize> = Vec::new(); // parties indexed from 0..n-1
    for i in 0..crs.n {
        if selector[i] {
            parties.push(i);
        } else {
            points.push(domain.roots_of_unity[i].clone());
        }
    }

    let b = interp_mostly_zero(&points);
    let b_evals = domain.fft(&b.coeffs, false).unwrap();

    debug_assert_eq!(
        b.degree() as usize,
        points.len(),
        "b.degree should be equal to points.len()"
    );
    debug_assert!(b.eval(&Scalar::zero()) == Scalar::one());

    // commit to b in g2
    let b_g2: V::Signature = crs.commit_sig(&b.coeffs);

    // q0 = (b-1)/x
    let q0_g1 = crs.compute_opening_proof(&b.coeffs, &Scalar::zero());

    // bhat = x^{t} * b
    // insert t 0s at the beginning of bhat.coeffs
    let mut bhat_coeffs = vec![Scalar::zero(); ct.t];
    bhat_coeffs.append(&mut b.coeffs.to_vec());
    let bhat = Poly::from_coeffs(NonEmptyVec::from_unchecked(bhat_coeffs));
    debug_assert_eq!(bhat.degree() as usize, crs.n);

    let bhat_g1: V::Public = crs.commit_public(&bhat.coeffs);

    let n_inv = Scalar::from_u64(crs.n as u64).inv();

    // compute the aggregate public key
    let mut bases: Vec<V::Public> = Vec::new();
    let mut scalars: Vec<Scalar> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.lag_pks[i].bls_pk);
        scalars.push(b_evals[i].clone());
    }
    let mut apk =
        <<V as Variant>::Public as Space<Scalar>>::msm(bases.as_slice(), scalars.as_slice(), 1);
    apk *= &n_inv;

    // compute sigma = (\sum B(omega^i)partial_decryptions[i])/(n) for i in parties
    let mut bases: Vec<V::Signature> = Vec::new();
    let mut scalars: Vec<Scalar> = Vec::new();
    for &i in &parties {
        bases.push(partial_decryptions[i].signature.into());
        scalars.push(b_evals[i].clone());
    }
    let mut sigma =
        <<V as Variant>::Signature as Space<Scalar>>::msm(bases.as_slice(), scalars.as_slice(), 1);
    sigma *= &n_inv;

    // compute Qx, Qhatx and Qz
    let mut bases: Vec<V::Public> = Vec::new();
    let mut scalars: Vec<Scalar> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.lag_pks[i].sk_li_x);
        scalars.push(b_evals[i].clone());
    }
    let qx =
        <<V as Variant>::Public as Space<Scalar>>::msm(bases.as_slice(), scalars.as_slice(), 1);

    let mut bases: Vec<V::Public> = Vec::new();
    let mut scalars: Vec<Scalar> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.lag_pks[i].sk_li_minus0);
        scalars.push(b_evals[i].clone());
    }
    let qhatx =
        <<V as Variant>::Public as Space<Scalar>>::msm(bases.as_slice(), scalars.as_slice(), 1);

    let mut bases: Vec<V::Public> = Vec::new();
    let mut scalars: Vec<Scalar> = Vec::new();
    for &i in &parties {
        bases.push(agg_key.agg_sk_li_lj_z[i]);
        scalars.push(b_evals[i].clone());
    }
    let qz =
        <<V as Variant>::Public as Space<Scalar>>::msm(bases.as_slice(), scalars.as_slice(), 1);

    // e(w1||sa1, sa2||w2)
    let minus1 = -Scalar::one();
    let w1 = [
        apk * &minus1,
        qz * &minus1,
        qx * &minus1,
        qhatx,
        bhat_g1 * &minus1,
        q0_g1 * &minus1,
    ];
    let w2 = [b_g2, sigma];

    let mut enc_key_lhs = w1.to_vec();
    enc_key_lhs.append(&mut ct.sa1.to_vec());

    let mut enc_key_rhs = ct.sa2.to_vec();
    enc_key_rhs.append(&mut w2.to_vec());

    // NOTE: Use multi pairings which are ~3x faster that doing individual pairings
    let enc_key = V::pairing(&w1[0], &ct.sa2[0])
        * &(V::pairing(&w1[1], &ct.sa2[1]))
        * &(V::pairing(&w1[2], &ct.sa2[2]))
        * &(V::pairing(&w1[3], &ct.sa2[3]))
        * &(V::pairing(&w1[4], &ct.sa2[4]))
        * &(V::pairing(&w1[5], &ct.sa2[5]))
        * &(V::pairing(&ct.sa1[0], &w2[0]))
        * &(V::pairing(&ct.sa1[1], &w2[1]));
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

    cipher.decrypt(&aes_nonce.into(), ct.ct.as_ref()).unwrap()
}

#[cfg(test)]
mod tests {
    use commonware_math::algebra::CryptoGroup;
    use commonware_math::algebra::Random;
    use rand::thread_rng;

    use crate::bls12381::{
        hints::{encryption::encrypt, setup::SecretKey},
        primitives::variant::{self, MinPk},
    };

    use super::*;

    #[test]
    fn test_decryption() {
        let mut rng = thread_rng();
        let n = 1 << 3; // actually n-1 total parties. one party is a dummy party that is always true
        let t: usize = 1;
        debug_assert!(t < n);

        let crs = CRS::new(n);

        let msg = b"Hello, world!";

        let sk = (0..n).map(|i| SecretKey::new(i)).collect::<Vec<_>>();

        let pk = sk
            .iter()
            .enumerate()
            .map(|(i, sk)| sk.get_lagrange_pk(i, &crs))
            .collect::<Vec<_>>();

        let (ak, ek) = AggregateKey::<MinPk>::new(pk, &crs);

        let gamma_g2 =
            <MinPk as variant::Variant>::Signature::generator() * &Scalar::random(&mut rng);
        let ct = encrypt::<MinPk>(&ek, t, &crs, gamma_g2, msg);

        // compute partial decryptions
        let mut partial_decryptions: Vec<PartialDecryption<MinPk>> = Vec::new();
        for i in 0..t {
            partial_decryptions.push(sk[i].partial_decryption(&ct));
        }
        for _ in t..n {
            partial_decryptions.push(PartialDecryption::zero());
        }

        // compute the decryption key
        let mut selector: Vec<bool> = Vec::new();
        for _ in 0..t {
            selector.push(true);
        }
        for _ in t..n {
            selector.push(false);
        }

        assert_eq!(
            agg_dec(&partial_decryptions, &ct, &selector, &ak, &crs),
            msg
        );
    }
}
