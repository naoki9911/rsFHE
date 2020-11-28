use crate::mulfft;
use crate::params;
use crate::tlwe;
use crate::utils;
use fftw::array::AlignedVec;
use fftw::types::*;
use rand::Rng;
use std::convert::TryInto;

#[derive(Debug, Copy, Clone)]
pub struct TRLWELv1 {
    pub a: [u32; params::trlwe_lv1::N],
    pub b: [u32; params::trlwe_lv1::N],
}

impl TRLWELv1 {
    pub fn new() -> TRLWELv1 {
        return TRLWELv1 {
            a: [0; params::trlwe_lv1::N],
            b: [0; params::trlwe_lv1::N],
        };
    }
}

pub fn trlweSymEncrypt(
    p: &Vec<f64>,
    alpha: f64,
    key: &Vec<u32>,
    twist: &AlignedVec<c64>,
) -> TRLWELv1 {
    let mut rng = rand::thread_rng();
    let mut trlwe = TRLWELv1::new();
    trlwe.a.iter_mut().for_each(|e| *e = rng.gen());

    let normal_distr = rand_distr::Normal::new(0.0, alpha).unwrap();
    let mut rng = rand::thread_rng();
    trlwe.b = utils::gussian_f64_vec(p, &normal_distr, &mut rng)
        .try_into()
        .unwrap();
    let a_i32 = trlwe.a.iter().map(|&e| e as i32).collect();
    let key_i32 = key.iter().map(|&e| e as i32).collect();
    let poly_res = mulfft::polynomial_mul(&a_i32, &key_i32, twist);

    for (bref, rval) in trlwe.b.iter_mut().zip(poly_res.iter()) {
        *bref = bref.wrapping_add(*rval);
    }

    return trlwe;
}

pub fn trlweSymDecrypt(trlwe: &TRLWELv1, key: &Vec<u32>, twist: &AlignedVec<c64>) -> Vec<u32> {
    let c_0_i32 = trlwe.a.iter().map(|&e| e as i32).collect();
    let key_i32 = key.iter().map(|&e| e as i32).collect();
    let poly_res = mulfft::polynomial_mul(&c_0_i32, &key_i32, twist);
    let mut res: Vec<u32> = Vec::new();
    for i in 0..trlwe.a.len() {
        let value = (trlwe.b[i].wrapping_sub(poly_res[i])) as i32;
        if value < 0 {
            res.push(0);
        } else {
            res.push(1);
        }
    }
    return res;
}

pub fn sample_extract_index(trlwe: &TRLWELv1, k: usize) -> tlwe::TLWELv1 {
    let mut res = tlwe::TLWELv1::new();

    const N: usize = params::trlwe_lv1::N;
    for i in 0..N {
        if i <= k {
            res.p[i] = trlwe.a[k - i];
        } else {
            res.p[i] = u32::MAX - trlwe.a[N + k - i];
        }
    }
    *res.b_mut() = trlwe.b[k];

    return res;
}

#[cfg(test)]
mod tests {
    use crate::mulfft;
    use crate::params;
    use crate::tlwe;
    use crate::trlwe;
    use rand::Rng;

    #[test]
    fn test_trlwe_enc_and_dec() {
        let mut rng = rand::thread_rng();

        // Generate 1024bits secret key
        let mut key: Vec<u32> = Vec::new();
        let mut key_dirty: Vec<u32> = Vec::new();
        const N: usize = params::trlwe_lv1::N;
        for i in 0..N {
            key.push((rng.gen::<u8>() % 2) as u32);
            key_dirty.push((rng.gen::<u8>() % 2) as u32);
        }

        let twist = mulfft::twist_gen(N);
        let mut correct = 0;
        let try_num = 500;

        for i in 0..try_num {
            let mut plain_text_enc: Vec<f64> = Vec::new();
            let mut plain_text: Vec<u32> = Vec::new();

            for j in 0..N {
                let sample: u32 = rng.gen::<u32>() % 2;
                let mut mu = 0.125;
                if sample == 0 {
                    mu = -0.125;
                }
                plain_text.push(sample);
                plain_text_enc.push(mu);
            }

            let c = trlwe::trlweSymEncrypt(&plain_text_enc, params::trlwe_lv1::ALPHA, &key, &twist);
            let dec = trlwe::trlweSymDecrypt(&c, &key, &twist);
            let dec_dirty = trlwe::trlweSymDecrypt(&c, &key_dirty, &twist);

            for j in 0..N {
                assert_eq!(plain_text[j], dec[j]);
                if plain_text[j] != dec_dirty[j] {
                    correct += 1;
                }
            }
        }

        let probability = correct as f64 / (try_num * N) as f64;
        assert!(probability - 0.50 < 0.1);
    }

    #[test]
    fn test_sample_extract_index() {
        let mut rng = rand::thread_rng();

        // Generate 1024bits secret key
        let mut key: Vec<u32> = Vec::new();
        let mut key_dirty: Vec<u32> = Vec::new();
        const N: usize = params::trlwe_lv1::N;
        for i in 0..N {
            key.push((rng.gen::<u8>() % 2) as u32);
            key_dirty.push((rng.gen::<u8>() % 2) as u32);
        }

        let twist = mulfft::twist_gen(N);
        let mut correct = 0;
        let try_num = 10;

        for i in 0..try_num {
            let mut plain_text_enc: Vec<f64> = Vec::new();
            let mut plain_text: Vec<u32> = Vec::new();

            for j in 0..N {
                let sample: u32 = rng.gen::<u32>() % 2;
                let mut mu = 0.125;
                if sample == 0 {
                    mu = -0.125;
                }
                plain_text.push(sample);
                plain_text_enc.push(mu);
            }

            let c = trlwe::trlweSymEncrypt(&plain_text_enc, params::trlwe_lv1::ALPHA, &key, &twist);

            for j in 0..N {
                let tlwe = trlwe::sample_extract_index(&c, j);
                let dec = tlwe::tlweLv1SymDecrypt(&tlwe, &key);
                let dec_dirty = tlwe::tlweLv1SymDecrypt(&tlwe, &key_dirty);
                assert_eq!(plain_text[j], dec);
                if plain_text[j] != dec_dirty {
                    correct += 1;
                }
            }
        }

        let probability = correct as f64 / (try_num * N) as f64;
        assert!(probability - 0.50 < 0.1);
    }
}
