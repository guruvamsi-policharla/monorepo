#[cfg(test)]
use crate::bls12381::primitives::group::Scalar;
#[cfg(test)]
use std::ops::{Add, Mul};

/// Simplified Discrete Fourier Transform, mainly used for testing
#[cfg(test)]
pub fn fft_fr_slow(
    ret: &mut [Scalar],
    data: &[Scalar],
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
    use commonware_math::algebra::Additive;

    use super::*;
    use crate::bls12381::{hints::fft_settings::Settings, primitives::group::Scalar};

    pub fn compare_sft_fft(
        fft_fr_slow: &dyn Fn(&mut [Scalar], &[Scalar], usize, &[Scalar], usize),
    ) {
        let size: usize = 12;

        let settings = Settings::new(size).unwrap();

        let data = (0..1 << size)
            .map(|i| Scalar::from_u64_arr(&[i as u64, 0, 0, 0]))
            .collect::<Vec<_>>();

        let mut out0 = vec![Scalar::zero(); settings.max_width];

        // Compare fast and slow FFT approach
        fft_fr_slow(&mut out0, &data, 1, &settings.expanded_roots_of_unity, 1);
        let out1 = settings.fft(&data, false).unwrap();

        for i in 0..settings.max_width {
            assert!(out0[i] == out1[i]);
        }
    }

    #[test]
    fn compare_sft_fft_() {
        compare_sft_fft(&fft_fr_slow);
    }

    #[test]
    fn fft_test() {
        let size: usize = 4;

        let settings = Settings::new(size).unwrap();

        println!("roots_of_unity: {:?}", settings.roots_of_unity);
        println!(
            "expanded_roots_of_unity: {:?}",
            settings.expanded_roots_of_unity
        );
        println!(
            "reverse_roots_of_unity: {:?}",
            settings.reverse_roots_of_unity
        );
        // let data = (0..1 << size)
        //     .map(|i| Scalar::from_u64_arr(&[i as u64, 0, 0, 0]))
        //     .collect::<Vec<_>>();

        // let mut out0 = vec![Scalar::zero(); settings.max_width];

        // // Compare fast and slow FFT approach
        // fft_fr_slow(&mut out0, &data, 1, &settings.expanded_roots_of_unity, 1);
        // let out1 = settings.fft(&data, false).unwrap();

        // println!("out1: {:?}", out1);
    }
}
