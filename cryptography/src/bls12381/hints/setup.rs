use commonware_math::algebra::{CryptoGroup, Space};
use commonware_math::{
    algebra::{Additive, Field, Random, Ring},
    poly::Poly,
};
use commonware_utils::vec::NonEmptyVec;
use rand::thread_rng;

use crate::bls12381::hints::types::Ciphertext;
use crate::bls12381::{
    hints::{crs::CRS, fft_settings::Settings, utils::lagrange_poly},
    primitives::{group::Scalar, variant::Variant},
};

#[derive(Clone)]
pub struct LagPolys {
    pub l: Vec<Poly<Scalar>>,
    pub l_minus0: Vec<Poly<Scalar>>,
    pub l_x: Vec<Poly<Scalar>>,
    pub li_lj_z: Vec<Vec<Poly<Scalar>>>,
    pub denom: Scalar,
}

impl LagPolys {
    // domain is the roots of unity of size n
    pub fn new(n: usize) -> Self {
        let domain = Settings::new(n.trailing_zeros() as usize).unwrap();
        // let domain = Radix2EvaluationDomain::<F>::new(n).unwrap();

        // compute polynomial L_i(X)
        let mut l = vec![Poly::zero(); n];
        for (i, ell) in l.iter_mut().enumerate().take(n) {
            *ell = lagrange_poly(n, i);
        }

        // compute polynomial (L_i(X) - L_i(0))*X
        let mut l_minus0 = vec![Poly::zero(); n];
        for i in 0..n {
            let mut li_minus0_coeffs = l[i].coeffs.clone();
            li_minus0_coeffs[0] = Scalar::zero();
            li_minus0_coeffs.insert(0, Scalar::zero());
            l_minus0[i] = Poly::from_coeffs(li_minus0_coeffs);
        }

        // compute polynomial (L_i(X) - L_i(0))/X
        let mut l_x = vec![Poly::zero(); n];
        for i in 0..n {
            l_x[i] = Poly::from_coeffs(NonEmptyVec::from_unchecked(
                l_minus0[i].coeffs[2..].to_vec(),
            ));
        }

        // compute polynomial L_i(X)*L_j(X)/Z(X) and (L_i(X)*L_i(X) - L_i(X))/Z(X)
        let mut li_lj_z = vec![vec![Poly::zero(); n]; n];
        for i in 0..n {
            for j in 0..n {
                li_lj_z[i][j] = if i == j {
                    (l[i].poly_mul(&l[i]).clone() - &l[i])
                        .divide_by_monomial(n, -Scalar::one())
                        .0
                    // divide_by_monomial(&(l[i].poly_mul(&l[i]).clone() - &l[i]), n, -Scalar::one()).0
                } else {
                    (l[i].poly_mul(&l[j]))
                        .divide_by_monomial(n, -Scalar::one())
                        .0
                    // (&l[i] * &l[j]).divide_by_vanishing_poly(domain).0
                };
            }
        }

        let mut denom = Scalar::one();
        for i in 1..n {
            denom *= &(Scalar::one() - &domain.roots_of_unity[i]);
        }

        // for i in 0..n {
        //     for j in 0..n {
        //         let monomial =
        //             DensePolynomial::from_coefficients_vec(vec![-domain.element(j), F::one()]);

        //         let computed = &l[i] / &monomial;
        //         assert_eq!(
        //             li_lj_z[i][j].evaluate(&F::zero()),
        //             computed.evaluate(&F::zero()) / (denom * domain.element(n - j))
        //         );
        //     }
        // }

        Self {
            l,
            l_minus0,
            l_x,
            li_lj_z,
            denom: denom.inv(),
        }
    }
}

#[derive(Clone)]
pub struct SecretKey {
    pub id: usize,
    sk: Scalar,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PartialDecryption<V: Variant> {
    /// Party id
    pub id: usize,
    /// Party commitment
    pub signature: V::Signature,
}

impl<V: Variant> PartialDecryption<V> {
    pub fn zero() -> Self {
        PartialDecryption {
            id: 0,
            signature: V::Signature::zero(),
        }
    }
}

/// Position oblivious public key -- slower to aggregate
#[derive(Clone)]
pub struct PublicKey<V: Variant> {
    pub bls_pk: V::Public,     //BLS pk
    pub hints: Vec<V::Public>, //hints
    // pub y: Vec<V::Public>, /* preprocessed toeplitz matrix. only for efficiency and can be
    //   * computed from hints */
    pub id: usize, // canonically assigned unique id in the system
}

/// Public key that can only be used in a fixed position -- faster to aggregate
#[derive(Clone)]
pub struct LagPublicKey<V: Variant> {
    pub id: usize,                  //id of the party
    pub position: usize,            //position in the aggregate key
    pub bls_pk: V::Public,          //BLS pk
    pub sk_li: V::Public,           //hint
    pub sk_li_minus0: V::Public,    //hint
    pub sk_li_lj_z: Vec<V::Public>, //hint
    pub sk_li_x: V::Public,         //hint
}

impl<V: Variant> LagPublicKey<V> {
    pub fn new(
        id: usize,
        position: usize,
        bls_pk: V::Public,
        sk_li: V::Public,
        sk_li_minus0: V::Public,
        sk_li_lj_z: Vec<V::Public>, //i = id
        sk_li_x: V::Public,
    ) -> Self {
        LagPublicKey {
            id,
            position,
            bls_pk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }
}

impl SecretKey {
    pub fn new(id: usize) -> Self {
        SecretKey {
            id,
            sk: Scalar::random(&mut thread_rng()),
        }
    }

    pub fn from_scalar(sk: Scalar, id: usize) -> Self {
        SecretKey { id, sk }
    }

    pub fn get_pk<V: Variant>(&self, crs: &CRS<V>) -> PublicKey<V> {
        let mut hints = vec![V::Public::zero(); crs.powers_of_g.len()];

        let bls_pk = V::Public::generator() * &self.sk;

        for (i, hint) in hints.iter_mut().enumerate().take(crs.powers_of_g.len()) {
            *hint = (crs.powers_of_g[i] * &self.sk).into();
        }

        // // compute y
        // let mut y = vec![V::Public::zero(); crs.y.len()];
        // for i in 0..crs.y.len() {
        //     y[i] = (crs.y[i] * &self.sk).into();
        // }

        PublicKey {
            id: self.id,
            bls_pk,
            hints,
            // y,
        }
    }

    pub fn get_lagrange_pk<V: Variant>(&self, position: usize, crs: &CRS<V>) -> LagPublicKey<V> {
        let mut sk_li_lj_z = vec![];

        let sk_li = crs.li[position] * &self.sk;

        let sk_li_minus0 = crs.li_minus0[position] * &self.sk;

        let sk_li_x = crs.li_x[position] * &self.sk;

        for j in 0..crs.n {
            sk_li_lj_z.push(crs.li_lj_z[position][j] * &self.sk);
        }

        LagPublicKey {
            id: self.id,
            position,
            bls_pk: V::Public::generator() * &self.sk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }

    pub fn partial_decryption<V: Variant>(&self, ct: &Ciphertext<V>) -> PartialDecryption<V> {
        PartialDecryption {
            id: self.id,
            signature: ct.gamma_g2 * &self.sk, // bls signature on gamma_g2
        }
    }
}

impl<V: Variant> PublicKey<V> {
    pub fn get_lag_public_key(
        &self,
        position: usize,
        crs: &CRS<V>,
        lag_polys: &LagPolys,
    ) -> LagPublicKey<V> {
        assert!(position < crs.n, "position out of bounds");

        let bls_pk = self.bls_pk;

        // compute sk_li
        let sk_li = <<V as Variant>::Public as Space<Scalar>>::msm(
            &self.hints[0..(lag_polys.l[position].degree() + 1) as usize],
            &*lag_polys.l[position].coeffs,
            1,
        );

        // compute sk_li_minus0
        let sk_li_minus0 = <<V as Variant>::Public as Space<Scalar>>::msm(
            &self.hints[0..(lag_polys.l_minus0[position].degree() + 1) as usize],
            &*lag_polys.l_minus0[position].coeffs,
            1,
        );

        // compute sk_li_x
        let sk_li_x = <<V as Variant>::Public as Space<Scalar>>::msm(
            &self.hints[0..(lag_polys.l_x[position].degree() + 1) as usize],
            &*lag_polys.l_x[position].coeffs,
            1,
        );

        // // compute sk*Li*Lj/Z = sk*Li/(X-omega^j)*(omega^j/denom) for all j in [n]\{i}
        // // for j = i: (Li^2 - Li)/Z = (Li - 1)/(X-omega^i)*(omega^i/denom)
        // // this is the same as computing KZG opening proofs at all points
        // // in the roots of unity domain for the polynomial Li(X), where the
        // // crs is {g^sk, g^{sk * tau}, g^{sk * tau^2}, ...}
        // // todo: move to https://eprint.iacr.org/2024/1279.pdf
        // let domain = Radix2EvaluationDomain::<E::ScalarField>::new(crs.n).unwrap();
        // let mut sk_li_lj_z = open_all_values::<E>(&self.y, &lag_polys.l[position].coeffs, &domain);

        // for (j, s) in sk_li_lj_z.iter_mut().enumerate().take(crs.n) {
        //     *s *= domain.element(j) * lag_polys.denom;
        // }

        // compute sk_li_lj_z
        let mut sk_li_lj_z = vec![V::Public::zero(); crs.n];

        for j in 0..crs.n {
            sk_li_lj_z[j] = <<V as Variant>::Public as Space<Scalar>>::msm(
                &self.hints[0..(lag_polys.li_lj_z[position][j].degree() + 1) as usize],
                &*lag_polys.li_lj_z[position][j].coeffs,
                1,
            );
        }

        // assert_eq!(sk_li_lj_z, my_sk_li_lj_z);

        LagPublicKey {
            id: self.id,
            position,
            bls_pk,
            sk_li,
            sk_li_minus0,
            sk_li_lj_z,
            sk_li_x,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bls12381::{hints::crs::CRS, primitives::variant::MinPk};

    #[test]
    fn test_setup() {
        let n = 1 << 4;
        let crs = CRS::<MinPk>::new(n);

        let mut sk: Vec<SecretKey> = Vec::new();
        let mut pk: Vec<LagPublicKey<MinPk>> = Vec::new();
        let mut lagrange_pk: Vec<LagPublicKey<MinPk>> = Vec::new();

        for i in 0..n {
            sk.push(SecretKey::new(i));
            pk.push(sk[i].get_lagrange_pk(i, &crs));
            lagrange_pk.push(sk[i].get_lagrange_pk(i, &crs));

            assert_eq!(pk[i].sk_li, lagrange_pk[i].sk_li);
            assert_eq!(pk[i].sk_li_minus0, lagrange_pk[i].sk_li_minus0);
            assert_eq!(pk[i].sk_li_x, lagrange_pk[i].sk_li_x);
            assert_eq!(pk[i].sk_li_lj_z, lagrange_pk[i].sk_li_lj_z);
        }

        // let _ak = AggregateKey::<MinPk>::new(pk, &crs);
    }

    #[test]
    fn test_setup_lag_setup() {
        let n = 1 << 4;
        let crs = CRS::<MinPk>::new(n);
        let lagpolys = LagPolys::new(n);

        let sk = SecretKey::new(0);
        let pk = sk.get_pk(&crs);
        let lag_pk = sk.get_lagrange_pk(0, &crs);

        let computed_lag_pk = pk.get_lag_public_key(0, &crs, &lagpolys);

        assert_eq!(computed_lag_pk.bls_pk, lag_pk.bls_pk);
        assert_eq!(computed_lag_pk.sk_li, lag_pk.sk_li);
        assert_eq!(computed_lag_pk.sk_li_minus0, lag_pk.sk_li_minus0);
        assert_eq!(computed_lag_pk.sk_li_x, lag_pk.sk_li_x);
        assert_eq!(computed_lag_pk.sk_li_lj_z, lag_pk.sk_li_lj_z);
    }
}
