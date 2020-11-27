use crate::mulfft;
use crate::trlwe;
use crate::tlwe;
use crate::utils;
use fftw::array::AlignedVec;
use fftw::types::*;

pub struct TRGSW {
    trlwe: Vec<trlwe::TRLWE>,
}

const fn N() -> usize {
    1024
}

const fn Nbit() -> usize {
    10
}

const fn BgBit() -> u32 {
    6
}

const fn Bg() -> u32 {
    1 << BgBit()
}

const fn l() -> usize {
    3
}

const fn alpha() -> f64 {
    2.98023223876953125e-08
}

pub fn trgswSymEncrypt(p: u32, alpha: f64, key: &Vec<u32>, twist: &AlignedVec<c64>) -> TRGSW {
    let mut p_f64: Vec<f64> = Vec::new();
    for i in 0..l() {
        let BgM: f64 = (Bg() as f64).powf(((1 + i) as f64) * -1.0);
        p_f64.push(BgM);
    }
    let p_torus = utils::f64_to_u32_torus(&p_f64);

    let mut plain_zero: Vec<f64> = Vec::new();
    for i in 0..N() {
        plain_zero.push(0.0);
    }

    let mut trgsw = TRGSW { trlwe: Vec::new() };

    for i in 0..l() * 2 {
        trgsw
            .trlwe
            .push(trlwe::trlweSymEncrypt(&plain_zero, alpha, &key, &twist));
    }

    for i in 0..l() {
        trgsw.trlwe[i].a[0] = trgsw.trlwe[i].a[0].wrapping_add(p * p_torus[i]);
        trgsw.trlwe[i + l()].b[0] = trgsw.trlwe[i + l()].b[0].wrapping_add(p * p_torus[i]);
    }
    return trgsw;
}

pub fn external_product(
    trgsw: &TRGSW,
    trlwe: &trlwe::TRLWE,
    twist: &AlignedVec<c64>,
) -> trlwe::TRLWE {
    let dec_a = decomposition(&trlwe.a);
    let dec_b = decomposition(&trlwe.b);
    let mut res: trlwe::TRLWE = trlwe::TRLWE {
        a: Vec::new(),
        b: Vec::new(),
    };

    for i in 0..N() {
        res.a.push(0);
        res.b.push(0);
    }

    for i in 0..l() {
        let tmp = mulfft::polynomial_mul_u32(&dec_a[i], &trgsw.trlwe[i].a, twist);
        for j in 0..N() {
            res.a[j] = res.a[j].wrapping_add(tmp[j]);
        }

        let tmp = mulfft::polynomial_mul_u32(&dec_b[i], &trgsw.trlwe[i + l()].a, twist);
        for j in 0..N() {
            res.a[j] = res.a[j].wrapping_add(tmp[j]);
        }

        let tmp = mulfft::polynomial_mul_u32(&dec_a[i], &trgsw.trlwe[i].b, twist);
        for j in 0..N() {
            res.b[j] = res.b[j].wrapping_add(tmp[j]);
        }

        let tmp = mulfft::polynomial_mul_u32(&dec_b[i], &trgsw.trlwe[i + l()].b, twist);
        for j in 0..N() {
            res.b[j] = res.b[j].wrapping_add(tmp[j]);
        }
    }

    return res;
}

pub fn decomposition(a: &Vec<u32>) -> Vec<Vec<u32>> {
    let offset = gen_offset();

    let mut a_tilda: Vec<u32> = Vec::new();
    let mut res: Vec<Vec<u32>> = Vec::new();

    for i in 0..N() {
        a_tilda.push(a[i as usize].wrapping_add(offset));
    }

    for i in 1..(l() + 1) as u32 {
        res.push(Vec::new());
        for j in 0..N() {
            let tmp =
                ((a_tilda[j as usize] >> (32 - BgBit() * i)) & (Bg() - 1)).wrapping_sub(Bg() / 2);
            res[(i - 1) as usize].push(tmp)
        }
    }

    return res;
}

// if cond == 0 then in1 else in2
pub fn cmux(
    in1: &trlwe::TRLWE,
    in2: &trlwe::TRLWE,
    cond: &TRGSW,
    twist: &AlignedVec<c64>,
) -> trlwe::TRLWE {
    let mut tmp = trlwe::TRLWE {
        a: Vec::new(),
        b: Vec::new(),
    };
    for i in 0..N() {
        tmp.a.push(in2.a[i].wrapping_sub(in1.a[i]));
        tmp.b.push(in2.b[i].wrapping_sub(in1.b[i]));
    }

    let tmp2 = external_product(cond, &tmp, twist);
    let mut res = trlwe::TRLWE {
        a: Vec::new(),
        b: Vec::new(),
    };
    for i in 0..N() {
        res.a.push(tmp2.a[i].wrapping_add(in1.a[i]));
        res.b.push(tmp2.b[i].wrapping_add(in1.b[i]));
    }

    return res;
}

pub fn gen_offset() -> u32 {
    let mut offset: u32 = 0;

    for i in 0..(l() as u32) {
        offset = offset.wrapping_add(Bg() / 2 * (1 << (32 - (i + 1) * BgBit())));
    }

    return offset;
}

pub fn blind_rotate(src : &tlwe::TLWELv0, testvec: &trlwe::TRLWE, bkey: &Vec<TRGSW>, twist:&AlignedVec<c64>) -> trlwe::TRLWE {
    let b_tilda = 2 * N() - (((src.b() as usize) + (1 << (31 - Nbit() -1))) >> (32-Nbit()-1));
    let mut res = trlwe::TRLWE{
        a:poly_mul_with_X_k(&testvec.a, b_tilda),
        b:poly_mul_with_X_k(&testvec.b, b_tilda),
    };

    let n = tlwe::n();
    for i in 0..n {
        let a_tilda = ((src.p[i as usize].wrapping_add((1<< (31 -Nbit()-1)))) >> (32-Nbit()-1)) as usize;
        let res2 = trlwe::TRLWE {
            a:poly_mul_with_X_k(&res.a, a_tilda),
            b:poly_mul_with_X_k(&res.b, a_tilda),
        };
        res = cmux(&res, &res2, &bkey[i as usize], twist);
    }

    return res;
}

pub fn poly_mul_with_X_k(a: &Vec<u32>, k:usize) -> Vec<u32> {
    let N = a.len();

    let mut res:Vec<u32> = Vec::new();
    for i in 0..N {
        res.push(0);
    }

    if k < N {
        for i in 0..(N - k) {
            res[i+k] = a[i];
        }
        for i in (N-k)..N {
            res[i+k-N] = u32::MAX - a[i];
        }
    } else {
        for i in 0..2*N - k {
            res[i+k-N] = u32::MAX - a[i];
        }
        for i in (2*N - k)..N {
            res[i-(2*N -k)] = a[i];
        }
    }

    return res;
}

pub fn generate_testvector() -> trlwe::TRLWE {
    let mut testvec = trlwe::TRLWE {
        a:Vec::new(),
        b:Vec::new(),
    };
    let b_vec:Vec<f64> = vec![0.125];
    let b_torus = utils::f64_to_u32_torus(&b_vec)[0];
    for i in 0..N() {
        testvec.a.push(0);
        testvec.b.push(b_torus);
    }

    return testvec;
}

#[cfg(test)]
mod tests {
    use crate::mulfft;
    use crate::trgsw::*;
    use crate::trlwe;
    use crate::tlwe;
    use crate::utils;
    use rand::Rng;

    #[test]
    fn test_decomposition() {
        let mut rng = rand::thread_rng();

        // Generate 1024bits secret key
        let mut key: Vec<u32> = Vec::new();
        for i in 0..N() {
            key.push((rng.gen::<u8>() % 2) as u32);
        }

        let twist = mulfft::twist_gen(N());
        let mut h: Vec<f64> = Vec::new();
        let try_num = 1000;

        for i in 1..l() + 1 {
            let tmp = (Bg() as f64).powf(-(i as f64));
            h.push(tmp);
        }

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

            let c = trlwe::trlweSymEncrypt(&plain_text_enc, alpha(), &key, &twist);
            let c_decomp_1 = decomposition(&c.a);
            let c_decomp_2 = decomposition(&c.b);
            let h_u32 = utils::f64_to_u32_torus(&h);
            let mut res = trlwe::TRLWE {
                a: Vec::new(),
                b: Vec::new(),
            };
            for j in 0..N() {
                let mut tmp0: u32 = 0;
                let mut tmp1: u32 = 0;
                for k in 0..l() {
                    tmp0 = tmp0.wrapping_add(c_decomp_1[k][j].wrapping_mul(h_u32[k]));
                    tmp1 = tmp1.wrapping_add(c_decomp_2[k][j].wrapping_mul(h_u32[k]));
                }
                res.a.push(tmp0);
                res.b.push(tmp1);
            }

            let dec = trlwe::trlweSymDecrypt(&res, &key, &twist);

            for j in 0..N() {
                assert_eq!(plain_text[j], dec[j]);
            }
        }
    }

    #[test]
    fn test_external_product() {
        let mut rng = rand::thread_rng();

        // Generate 1024bits secret key
        let mut key: Vec<u32> = Vec::new();
        for i in 0..N() {
            key.push((rng.gen::<u8>() % 2) as u32);
        }

        let twist = mulfft::twist_gen(N());
        let try_num = 100;

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

            let c = trlwe::trlweSymEncrypt(&plain_text_enc, alpha(), &key, &twist);
            let p = trlwe::trlweSymDecrypt(&c, &key, &twist);
            let trgsw_true = trgswSymEncrypt(1, alpha(), &key, &twist);
            let ext_c = external_product(&trgsw_true, &c, &twist);
            let dec = trlwe::trlweSymDecrypt(&ext_c, &key, &twist);

            for j in 0..N() {
                assert_eq!(plain_text[j], p[j]);
            }
            for j in 0..N() {
                assert_eq!(plain_text[j], dec[j]);
            }
        }
    }

    #[test]
    fn test_cmux() {
        let mut rng = rand::thread_rng();
        let mut key: Vec<u32> = Vec::new();
        for i in 0..N() {
            key.push((rng.gen::<u8>() % 2) as u32);
        }

        let twist = mulfft::twist_gen(N());
        let try_num = 100;
        for i in 0..try_num {
            let mut plain_text_enc_1: Vec<f64> = Vec::new();
            let mut plain_text_enc_2: Vec<f64> = Vec::new();
            let mut plain_text_1: Vec<u32> = Vec::new();
            let mut plain_text_2: Vec<u32> = Vec::new();

            for j in 0..N() {
                let sample: u32 = rng.gen::<u32>() % 2;
                let mut mu = 0.125;
                if sample == 0 {
                    mu = -0.125;
                }
                plain_text_1.push(sample);
                plain_text_enc_1.push(mu);
            }
            for j in 0..N() {
                let sample: u32 = rng.gen::<u32>() % 2;
                let mut mu = 0.125;
                if sample == 0 {
                    mu = -0.125;
                }
                plain_text_2.push(sample);
                plain_text_enc_2.push(mu);
            }
            let c1 = trlwe::trlweSymEncrypt(&plain_text_enc_1, alpha(), &key, &twist);
            let c2 = trlwe::trlweSymEncrypt(&plain_text_enc_2, alpha(), &key, &twist);
            let trgsw_true = trgswSymEncrypt(1, alpha(), &key, &twist);
            let trgsw_false = trgswSymEncrypt(0, alpha(), &key, &twist);
            let enc_1 = cmux(&c1, &c2, &trgsw_false, &twist);
            let enc_2 = cmux(&c1, &c2, &trgsw_true, &twist);
            let dec_1 = trlwe::trlweSymDecrypt(&enc_1, &key, &twist);
            let dec_2 = trlwe::trlweSymDecrypt(&enc_2, &key, &twist);
            for j in 0..N() {
                assert_eq!(plain_text_1[j], dec_1[j]);
                assert_eq!(plain_text_2[j], dec_2[j]);
            }
        }
    }

    #[test]
    fn test_blind_rotate() {
        let mut rng = rand::thread_rng();
        let twist = mulfft::twist_gen(N());
        let mut keyLv0: Vec<u32> = Vec::new();
        let mut keyLv1: Vec<u32> = Vec::new();
        for i in 0..tlwe::n() {
            keyLv0.push((rng.gen::<u8>() % 2) as u32);
        }
        for i in 0..N() {
            keyLv1.push((rng.gen::<u8>() % 2) as u32);
        }

        let mut bKey : Vec<TRGSW> = Vec::new();
        for i in 0..keyLv0.len() {
            bKey.push(trgswSymEncrypt(keyLv0[i], trlwe::alpha(), &keyLv1, &twist));
        }

        let try_num = 10;
        let test_vec = generate_testvector();
        for i in 0..try_num {
            let plain_text = rng.gen::<u32>() % 2;
            let mut mu = 0.125;
            if plain_text == 0 {
                mu = -0.125;
            }

            let tlwe = tlwe::tlweSymEncrypt(mu, tlwe::alpha(), &keyLv0);
            let trlwe = blind_rotate(&tlwe, &test_vec, &bKey, &twist);
            let tlwe_lv1 = trlwe::sample_extract_index(&trlwe, 0);
            let dec = tlwe::tlweLv1SymDecrypt(&tlwe_lv1, &keyLv1);
            assert_eq!(plain_text, dec);
        }
    }
}
