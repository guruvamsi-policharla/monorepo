#[cfg(test)]
use std::ops::{Add, Mul};

#[cfg(test)]
use crate::bls12381::primitives::group::{Scalar, G1};

// Used for testing
#[cfg(test)]
pub fn fft_g1_slow(
    ret: &mut [G1],
    data: &[G1],
    stride: usize,
    roots: &[Scalar],
    roots_stride: usize,
) {
    for i in 0..data.len() {
        // Evaluate first member at 1
        ret[i] = data[0].clone().mul(&roots[0]);

        // Evaluate the rest of members using a step of (i * J) % data.len() over the roots
        // This distributes the roots over correct x^n members and saves on multiplication
        for j in 1..data.len() {
            let v = data[j * stride]
                .clone()
                .mul(&roots[((i * j) % data.len()) * roots_stride]);
            ret[i] = ret[i].clone().add(&v);
        }
    }
}

#[cfg(test)]
mod tests {
    use commonware_math::algebra::{Additive, CryptoGroup};
    use std::ops::Mul;

    use super::*;
    use crate::bls12381::hints::fft_settings::Settings;
    use crate::bls12381::primitives::group::{Scalar, G1};

    pub fn compare_fft_g1(fft_g1_slow: &dyn Fn(&mut [G1], &[G1], usize, &[Scalar], usize)) {
        let size: usize = 6;
        let settings = Settings::new(size).unwrap();

        let data = (0..settings.max_width)
            .map(|i| G1::generator().mul(&Scalar::from_u64((i as u64) + 1)))
            .collect::<Vec<_>>();

        let mut out_slow = vec![G1::zero(); settings.max_width];

        fft_g1_slow(
            &mut out_slow,
            &data,
            1,
            &settings.expanded_roots_of_unity,
            1,
        );

        let out_fast = settings.fft_g1(&data, false).unwrap();

        for i in 0..settings.max_width {
            assert_eq!(out_slow[i], out_fast[i]);
        }
    }

    #[test]
    fn compare_fft_g1_() {
        compare_fft_g1(&fft_g1_slow);
    }
}
