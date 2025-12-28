use commonware_math::algebra::Additive;
use commonware_math::poly::Poly;

use crate::bls12381::primitives::group::Scalar;

/// Divide the polynomial `poly` by the monomial term `(X^degree + constant)`.
/// Returns `(quotient, remainder)` such that `poly = quotient * (X^degree + constant) + remainder`.
/// This is a generalized implementation of synthetic division and runs in O(deg(poly)) time.
pub fn divide_by_monomial(
    poly: &Poly<Scalar>,
    degree: usize,
    constant: &Scalar,
) -> (Poly<Scalar>, Poly<Scalar>) {
    let poly_degree = poly.degree() as usize;
    // If degree of poly < divisor degree, quotient is 0, remainder is poly.
    if poly_degree < degree {
        return (Poly::from_coeffs(vec![Scalar::zero()]), poly.clone());
    }

    // Quotient will have degree (n - 1) - degree = n - 1 - degree.
    // Size = n - degree.
    let mut quotient = vec![Scalar::zero(); poly_degree + 1 - degree];
    let mut remainder = poly.clone();

    // Synthetic division from high degree to low
    for i in (degree..poly_degree + 1).rev() {
        let q = remainder[i];
        quotient[i - degree] = q;

        let term = &q * constant;
        remainder[i - degree] -= &term;
    }

    // Truncate remainder to degree < degree (size degree).
    remainder.truncate(degree);
    // If remainder is empty (degree=0 case?), make it zero.
    if remainder.is_empty() {
        remainder.push(Scalar::zero());
    }

    (Poly::from_coeffs(quotient), Poly::from_coeffs(remainder))
}

// 1 at omega^i and 0 elsewhere on domain {omega^i}_{i \in [n]}
// Computed by dividing X^n - 1 by X - omega^i
pub fn lagrange_poly(n: usize, i: usize) -> Poly<Scalar> {
    debug_assert!(i < n);

    let mut vanishing_poly = Poly::from_coeffs(vec![Scalar::zero(); n + 1]);
    vanishing_poly[n] = Scalar::one();

    // divide by monomial (X - omega^i)
}

/// interpolates a polynomial which is zero on points and 1 at the point 0
/// todo: use faster interpolation
pub fn interp_mostly_zero<F: Field>(points: &Vec<F>) -> DensePolynomial<F> {
    if points.is_empty() {
        // threshold=n
        return DensePolynomial::from_coefficients_vec(vec![F::one()]);
    }

    let mut interp = DensePolynomial::from_coefficients_vec(vec![F::one()]);
    for &point in points {
        interp = interp.naive_mul(&DensePolynomial::from_coefficients_vec(vec![
            -point,
            F::one(),
        ]));
    }

    let scale = interp.evaluate(&F::zero());
    interp = &interp * (F::one() / scale);

    interp
}

/// Computes all the openings of a KZG commitment in O(n log n) time
/// See https://github.com/khovratovich/Kate/blob/master/Kate_amortized.pdf
/// eprint version has a bug and hasn't been updated
pub fn open_all_values<E: Pairing>(
    y: &[E::G1Affine],
    f: &[E::ScalarField],
    domain: &Radix2EvaluationDomain<E::ScalarField>,
) -> Vec<E::G1> {
    let top_domain = Radix2EvaluationDomain::<E::ScalarField>::new(2 * domain.size()).unwrap();

    // use FK22 to get all the KZG proofs in O(nlog n) time =======================
    // f = {f0 ,f1, ..., fd}
    // v = {(d 0s), f1, ..., fd}
    let mut v = vec![E::ScalarField::zero(); domain.size() + 1];
    v.append(&mut f[1..f.len()].to_vec());

    debug_assert_eq!(v.len(), 2 * domain.size());
    let v = top_domain.fft(&v);

    // h = y \odot v
    let mut h = vec![E::G1::zero(); 2 * domain.size()];
    for i in 0..2 * domain.size() {
        h[i] = y[i] * (v[i]);
    }

    // inverse fft on h
    let mut h = top_domain.ifft(&h);

    h.truncate(domain.size());

    // fft on h to get KZG proofs
    domain.fft(&h)
}

/// interpolates a polynomial where evaluations on points are zero and the polynomial evaluates to 1
/// at the point 1 but relies on the number of points being a power of 2
/// currently not used as this portion is not a bottleneck
// pub fn compute_vanishing_poly(points: &Vec<ScalarField>) -> DensePolynomial {
//     let mut monomials = Vec::new();
//     for i in 0..points.len() {
//         monomials.push(DensePolynomial::from_coeffs(
//             HostSlice::from_slice(&vec![ScalarField::zero() - points[i], ScalarField::one()]),
//             2,
//         ));
//     }

//     // assert that points.len() is a power of 2
//     assert_eq!(
//         points.len().count_ones(),
//         1,
//         "Implementation demands that n-t is a power of 2. Currently: {}",
//         points.len()
//     );

//     let mut chunk_size = points.len() / 2;
//     while chunk_size > 0 {
//         for i in 0..chunk_size {
//             monomials[i] = &monomials[i] * &monomials[i + chunk_size];
//         }
//         chunk_size = chunk_size / 2;
//     }

//     let scale = monomials[0].eval(&ScalarField::one());
//     let res = &monomials[0] * &scale.inv();

//     res
// }

#[cfg(test)]
mod tests {
    use ark_bls12_381::Bls12_381;
    use ark_ec::{bls12::Bls12, pairing::Pairing, VariableBaseMSM};
    use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
    use ark_std::{UniformRand, Zero};

    use crate::crs::CRS;

    use super::*;
    type Fr = <Bls12<ark_bls12_381::Config> as Pairing>::ScalarField;
    type G1 = <Bls12<ark_bls12_381::Config> as Pairing>::G1;
    type E = Bls12_381;

    #[test]
    fn open_all_test() {
        let mut rng = ark_std::test_rng();

        let n = 1 << 8;
        let domain = Radix2EvaluationDomain::<Fr>::new(n).unwrap();
        let crs = CRS::<E>::new(n, &mut ark_std::test_rng());

        let mut f: Vec<ark_ff::Fp<ark_ff::MontBackend<ark_bls12_381::FrConfig, 4>, 4>> =
            vec![Fr::zero(); n];
        for i in 0..n {
            f[i] = Fr::rand(&mut rng);
        }

        let com = G1::msm(&crs.powers_of_g[0..f.len()], &f).unwrap();

        let timer = std::time::Instant::now();
        let pi = open_all_values::<E>(&crs.y, &f, &domain);
        println!("open_all_values took {:?}", timer.elapsed());

        // verify the kzg proof
        let g = crs.powers_of_g[0];
        let h = crs.powers_of_h[0];

        let fpoly = DensePolynomial::from_coefficients_vec(f.clone());
        for i in 0..n {
            let lhs = E::pairing(com - (g * fpoly.evaluate(&domain.element(i))), h);
            let rhs = E::pairing(pi[i], crs.powers_of_h[1] - (h * domain.element(i)));
            assert_eq!(lhs, rhs);
        }
    }
}
