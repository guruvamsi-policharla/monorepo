use commonware_math::{
    algebra::{Additive, CryptoGroup, Field, Ring, Space},
    poly::Poly,
};
use commonware_utils::vec::NonEmptyVec;

use crate::bls12381::{
    hints::{fft_settings::Settings, utils::lagrange_poly},
    primitives::{group::Scalar, variant::Variant},
};

#[derive(Clone, Debug)]
pub struct CRS<V: Variant> {
    pub n: usize, // maximum number of parties in a committee
    pub powers_of_g: Vec<V::Public>,
    pub powers_of_h: Vec<V::Signature>,

    // preprocessed lagrange polynomials
    pub li: Vec<V::Public>,
    pub li_minus0: Vec<V::Public>,
    pub li_x: Vec<V::Public>,
    pub li_lj_z: Vec<Vec<V::Public>>,

    // preprocessed lagrange polynomials in g2 (only needed for verifying hints)
    pub li_g2: Vec<V::Signature>,
    pub li_minus0_g2: Vec<V::Signature>,
    pub li_x_g2: Vec<V::Signature>,
    pub li_lj_z_g2: Vec<Vec<V::Signature>>,

    // preprocessed Toeplitz matrix
    pub y: Vec<V::Public>,
}

impl<V: Variant> CRS<V> {
    pub fn new(n: usize) -> Self {
        // let tau = Scalar::random(&mut thread_rng());
        let tau = Scalar::from_u64(2);
        Self::deterministic_new(n, tau)
    }

    pub fn deterministic_new(n: usize, tau: Scalar) -> Self {
        let mut powers_of_tau = vec![Scalar::zero(); n + 1];
        let mut powers_of_g = vec![V::Public::generator(); n + 1];
        let mut powers_of_h = vec![V::Signature::generator(); n + 1];

        let mut cur = Scalar::one();
        for i in 0..=n {
            powers_of_tau[i] = cur.clone();
            powers_of_g[i] = V::Public::generator() * &cur;
            powers_of_h[i] = V::Signature::generator() * &cur;
            cur *= &tau;
        }

        // lagrange powers
        let mut li_evals: Vec<Scalar> = vec![Scalar::zero(); n];
        let mut li_evals_minus0: Vec<Scalar> = vec![Scalar::zero(); n];
        let mut li_evals_x: Vec<Scalar> = vec![Scalar::zero(); n];

        let tau2_inv: Scalar = (tau.clone() * &tau).inv();
        for i in 0..n {
            let li = lagrange_poly(n, i);
            li_evals[i] = li.eval(&tau);

            li_evals_minus0[i] = (li_evals[i].clone() - &li.coeffs[0]) * &tau;

            li_evals_x[i] = li_evals_minus0[i].clone() * &tau2_inv;
        }

        let tau_n = (0..n).fold(Scalar::one(), |acc, _| acc * &tau);

        let z_eval = tau_n - &Scalar::one();
        let z_eval_inv = z_eval.inv();

        let mut li = vec![V::Public::zero(); n];
        let mut li_g2 = vec![V::Signature::zero(); n];

        for i in 0..n {
            li[i] = V::Public::generator() * &li_evals[i];
            li_g2[i] = V::Signature::generator() * &li_evals[i];
        }

        let mut li_minus0 = vec![V::Public::zero(); n];
        let mut li_minus0_g2 = vec![V::Signature::zero(); n];

        for i in 0..n {
            li_minus0[i] = V::Public::generator() * &li_evals_minus0[i];
            li_minus0_g2[i] = V::Signature::generator() * &li_evals_minus0[i];
        }

        let mut li_x = vec![V::Public::zero(); n];
        let mut li_x_g2 = vec![V::Signature::zero(); n];

        for i in 0..n {
            li_x[i] = V::Public::generator() * &li_evals_x[i];
            li_x_g2[i] = V::Signature::generator() * &li_evals_x[i];
        }

        let mut li_lj_z = vec![vec![V::Public::zero(); n]; n];
        let mut li_lj_z_g2 = vec![vec![V::Signature::zero(); n]; n];

        for i in 0..n {
            for j in 0..n {
                li_lj_z[i][j] = if i == j {
                    V::Public::generator()
                        * &((li_evals[i].clone() * &li_evals[i] - &li_evals[i]) * &z_eval_inv)
                } else {
                    V::Public::generator() * &(li_evals[i].clone() * &li_evals[j] * &z_eval_inv)
                };

                li_lj_z_g2[i][j] = if i == j {
                    V::Signature::generator()
                        * &((li_evals[i].clone() * &li_evals[i] - &li_evals[i]) * &z_eval_inv)
                } else {
                    V::Signature::generator() * &(li_evals[i].clone() * &li_evals[j] * &z_eval_inv)
                };
            }
        }

        // Compute the Toeplitz matrix preprocessing
        // ==================================================
        let mut top_tau = powers_of_tau.clone();
        top_tau.truncate(n);
        top_tau.reverse();
        top_tau.resize(2 * n, Scalar::zero());

        let top_domain = Settings::new((2 * n).trailing_zeros() as usize).unwrap();
        let top_tau = top_domain.fft(&top_tau, false).unwrap();

        // Compute powers of top_tau
        let y = top_tau
            .iter()
            .map(|x| V::Public::generator() * x)
            .collect::<Vec<_>>();

        Self {
            n,
            powers_of_g,
            powers_of_h,

            li,
            li_minus0,
            li_x,
            li_lj_z,

            li_g2,
            li_minus0_g2,
            li_x_g2,
            li_lj_z_g2,

            y,
        }
    }

    pub fn commit_public(&self, coeffs: &[Scalar]) -> V::Public {
        assert!(
            coeffs.len() <= self.powers_of_g.len(),
            "Too many coefficients for the given powers of tau"
        );

        <<V as Variant>::Public as Space<Scalar>>::msm(&self.powers_of_g[..coeffs.len()], coeffs, 1)
    }

    pub fn commit_sig(&self, coeffs: &[Scalar]) -> V::Signature {
        assert!(
            coeffs.len() <= self.powers_of_h.len(),
            "Too many coefficients for the given powers of tau"
        );

        <<V as Variant>::Signature as Space<Scalar>>::msm(
            &self.powers_of_h[..coeffs.len()],
            coeffs,
            1,
        )
    }

    pub fn compute_opening_proof(&self, coeffs: &[Scalar], point: &Scalar) -> V::Public {
        let polynomial = Poly::from_coeffs(NonEmptyVec::from_unchecked(coeffs.to_vec()));
        let witness_polynomial = polynomial.divide_by_monomial(1, -point.clone()).0;

        self.commit_public(&witness_polynomial.coeffs)
    }
}

#[cfg(test)]
mod tests {
    use commonware_math::{
        algebra::{Additive, CryptoGroup, Field, Random, Ring},
        poly::Poly,
    };
    use commonware_utils::vec::NonEmptyVec;
    use rand::thread_rng;

    use crate::bls12381::{
        hints::{
            crs::{self, CRS},
            fft_settings::Settings,
        },
        primitives::{
            group::Scalar,
            variant::{self, MinPk, Variant},
        },
    };

    #[test]
    fn test_kzg() {
        let n = 1 << 3;
        let crs: CRS<MinPk> = crs::CRS::new(n);

        // sample n random coeffs
        let coeffs = (0..n)
            .map(|_| Scalar::random(&mut thread_rng()))
            .collect::<Vec<_>>();
        let com = crs.commit_public(&coeffs);

        let point = Scalar::random(&mut thread_rng());
        let eval = Poly::from_coeffs(NonEmptyVec::from_unchecked(coeffs.clone())).eval(&point);

        let pi = crs.compute_opening_proof(&coeffs, &point);

        let lhs = MinPk::pairing(
            &(com + &(<variant::MinPk as variant::Variant>::Public::generator() * &(-eval))),
            &<variant::MinPk as variant::Variant>::Signature::generator(),
        );
        let rhs = MinPk::pairing(
            &pi,
            &(crs.powers_of_h[1]
                - &(<variant::MinPk as variant::Variant>::Signature::generator() * &point)),
        );
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn test_sumcheck() {
        // A(X).B(X) = \sum_i A(i).B(i) + X * Q_x(X) + Z(X) * Q_Z(X)
        let n: usize = 1 << 5;
        let domain = Settings::new((n).trailing_zeros() as usize).unwrap();

        // sample n random evals
        let a_evals = (0..n)
            .map(|_| Scalar::random(&mut thread_rng()))
            .collect::<Vec<_>>();
        let b_evals = (0..n)
            .map(|_| Scalar::random(&mut thread_rng()))
            .collect::<Vec<_>>();
        let mut s = Scalar::zero();
        for i in 0..n {
            s += &(a_evals[i].clone() * &b_evals[i]);
        }

        let a_coeffs = domain.fft(&a_evals, true).unwrap();
        let b_coeffs = domain.fft(&b_evals, true).unwrap();

        let a_poly = Poly::from_coeffs(NonEmptyVec::from_unchecked(a_coeffs));
        let b_poly = Poly::from_coeffs(NonEmptyVec::from_unchecked(b_coeffs));

        let c_poly = a_poly.poly_mul(&b_poly);

        println!("a_poly deg: {}", a_poly.degree());
        println!("b_poly deg: {}", b_poly.degree());
        println!("c_poly deg: {}", c_poly.degree());

        let (qz, rem) = c_poly.divide_by_monomial(n, -Scalar::one());

        println!("qz deg: {}", qz.degree());
        println!("rem deg: {}", rem.degree());

        assert_eq!(
            s * &Scalar::from_u64(n as u64).inv(),
            rem.eval(&Scalar::zero())
        );
    }
}
