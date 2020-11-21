use crate::mulfft;
use crate::utils;
use fftw::array::AlignedVec;
use fftw::types::*;
use rand::Rng;

pub struct TRLWE {
    pub a: Vec<u32>,
    pub b: Vec<u32>,
}

pub struct TLWELv1 {
    pub p: Vec<u32>,
}

const fn N() -> usize {
    1024
}

const fn alpha() -> f64 {
    2.98023223876953125e-08
}

pub fn trlweSymEncrypt(p: &Vec<f64>, alpha: f64, key: &Vec<u32>, twist: &AlignedVec<c64>) -> TRLWE {
    let mut rng = rand::thread_rng();
    let mut trlwe: TRLWE = TRLWE {
        a: Vec::new(),
        b: Vec::new(),
    };
    for i in 0..key.len() {
        trlwe.a.push(rng.gen());
    }
    trlwe.b = utils::gussian_32bit(&utils::f64_to_u32_torus(p), alpha, key.len());
    let a_i32 = trlwe.a.iter().map(|&e| e as i32).collect();
    let key_i32 = key.iter().map(|&e| e as i32).collect();
    let poly_res = mulfft::polynomial_mul(&a_i32, &key_i32, twist);
    for i in 0..trlwe.b.len() {
        trlwe.b[i] = trlwe.b[i].wrapping_add(poly_res[i]);
    }

    return trlwe;
}

pub fn trlweSymDecrypt(trlwe: &TRLWE, key: &Vec<u32>, twist: &AlignedVec<c64>) -> Vec<u32> {
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

pub fn sample_extract_index(trlwe: &TRLWE, k: usize) -> TLWELv1 {
    let mut res = TLWELv1 {
        p : Vec::new()
    };

    for i in 0..N() {
        if i <= k {
            res.p.push(trlwe.a[k-i]);
        }else{
            res.p.push(u32::MAX-trlwe.a[N()+k-i]);
        }
    }
    res.p.push(trlwe.b[k]);

    return res;
}

#[cfg(test)]
mod tests {
    use crate::mulfft;
    use crate::trlwe::*;
    use crate::tlwe::*;

    #[test]
    fn test_trlwe_enc_and_dec() {
        let mut rng = rand::thread_rng();

        // Generate 1024bits secret key
        let mut key: Vec<u32> = Vec::new();
        let mut key_dirty: Vec<u32> = Vec::new();
        for i in 0..N() {
            key.push((rng.gen::<u8>() % 2) as u32);
            key_dirty.push((rng.gen::<u8>() % 2) as u32);
        }

        let twist = mulfft::twist_gen(N());
        let mut correct = 0;
        let try_num = 500;

        for i in 0..try_num {
            let mut plain_text_enc: Vec<f64> = Vec::new();
            let mut plain_text: Vec<u32> = Vec::new();

            for j in 0..N() {
                let sample: u32 = rng.gen::<u32>() % 2;
                let mut mu = 0.125;
                if sample == 0 {
                    mu = -0.125;
                }
                plain_text.push(sample);
                plain_text_enc.push(mu);
            }

            let c = trlweSymEncrypt(&plain_text_enc, alpha(), &key, &twist);
            let dec = trlweSymDecrypt(&c, &key, &twist);
            let dec_dirty = trlweSymDecrypt(&c, &key_dirty, &twist);

            for j in 0..N() {
                assert_eq!(plain_text[j], dec[j]);
                if plain_text[j] != dec_dirty[j] {
                    correct += 1;
                }
            }
        }

        let probability = correct as f64 / (try_num * N()) as f64;
        assert!(probability - 0.50 < 0.1);
    }

    #[test]
    fn test_sample_extract_index() {
        let mut rng = rand::thread_rng();

        // Generate 1024bits secret key
        let mut key: Vec<u32> = Vec::new();
        let mut key_dirty: Vec<u32> = Vec::new();
        for i in 0..N() {
            key.push((rng.gen::<u8>() % 2) as u32);
            key_dirty.push((rng.gen::<u8>() % 2) as u32);
        }

        let twist = mulfft::twist_gen(N());
        let mut correct = 0;
        let try_num = 10;

        for i in 0..try_num {
            let mut plain_text_enc: Vec<f64> = Vec::new();
            let mut plain_text: Vec<u32> = Vec::new();

            for j in 0..N() {
                let sample: u32 = rng.gen::<u32>() % 2;
                let mut mu = 0.125;
                if sample == 0 {
                    mu = -0.125;
                }
                plain_text.push(sample);
                plain_text_enc.push(mu);
            }

            let c = trlweSymEncrypt(&plain_text_enc, alpha(), &key, &twist);

            for j in 0..N() {
                let tlwe = sample_extract_index(&c, j);
                let dec = tlweLv1SymDecrypt(&tlwe, &key);
                let dec_dirty = tlweLv1SymDecrypt(&tlwe, &key_dirty);
                assert_eq!(plain_text[j], dec);
                if plain_text[j] != dec_dirty {
                    correct += 1;
                }
            }
        }

        let probability = correct as f64 / (try_num * N()) as f64;
        assert!(probability - 0.50 < 0.1);
    }
}
