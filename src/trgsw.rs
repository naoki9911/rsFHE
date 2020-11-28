use crate::mulfft;
use crate::params;
use crate::tlwe;
use crate::trlwe;
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

const fn basebit() -> usize {
    2
}

const fn iks_t() -> usize {
    8
}

pub fn trgswSymEncrypt(p: u32, alpha: f64, key: &Vec<u32>, twist: &AlignedVec<c64>) -> TRGSW {
    let mut p_f64: Vec<f64> = Vec::new();
    for i in 0..l() {
        let BgM: f64 = (Bg() as f64).powf(((1 + i) as f64) * -1.0);
        p_f64.push(BgM);
    }
    let p_torus = utils::f64_to_torus_vec(&p_f64);

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

pub fn blind_rotate(
    src: &tlwe::TLWELv0,
    testvec: &trlwe::TRLWE,
    bkey: &Vec<TRGSW>,
    twist: &AlignedVec<c64>,
) -> trlwe::TRLWE {
    let b_tilda = 2 * N() - (((src.b() as usize) + (1 << (31 - Nbit() - 1))) >> (32 - Nbit() - 1));
    let mut res = trlwe::TRLWE {
        a: poly_mul_with_X_k(&testvec.a, b_tilda),
        b: poly_mul_with_X_k(&testvec.b, b_tilda),
    };

    for i in 0..params::tlwe_lv0::N {
        let a_tilda = ((src.p[i as usize].wrapping_add((1 << (31 - Nbit() - 1))))
            >> (32 - Nbit() - 1)) as usize;
        let res2 = trlwe::TRLWE {
            a: poly_mul_with_X_k(&res.a, a_tilda),
            b: poly_mul_with_X_k(&res.b, a_tilda),
        };
        res = cmux(&res, &res2, &bkey[i as usize], twist);
    }

    return res;
}

pub fn poly_mul_with_X_k(a: &Vec<u32>, k: usize) -> Vec<u32> {
    let N = a.len();

    let mut res: Vec<u32> = Vec::new();
    for i in 0..N {
        res.push(0);
    }

    if k < N {
        for i in 0..(N - k) {
            res[i + k] = a[i];
        }
        for i in (N - k)..N {
            res[i + k - N] = u32::MAX - a[i];
        }
    } else {
        for i in 0..2 * N - k {
            res[i + k - N] = u32::MAX - a[i];
        }
        for i in (2 * N - k)..N {
            res[i - (2 * N - k)] = a[i];
        }
    }

    return res;
}

pub fn generate_testvector() -> trlwe::TRLWE {
    let mut testvec = trlwe::TRLWE {
        a: Vec::new(),
        b: Vec::new(),
    };
    let b_torus = utils::f64_to_torus(0.125);
    for i in 0..N() {
        testvec.a.push(0);
        testvec.b.push(b_torus);
    }

    return testvec;
}

pub fn identity_key_switching(
    src: &tlwe::TLWELv1,
    ks: &Vec<Vec<Vec<tlwe::TLWELv0>>>,
) -> tlwe::TLWELv0 {
    let mut res = tlwe::TLWELv0::new();

    res.p[params::tlwe_lv0::N] = src.p[src.p.len() - 1];

    let prec_offset = 1 << (32 - (1 + basebit() * iks_t()));

    for i in 0..N() {
        let ks_i = &ks[i];
        let a_bar = src.p[i].wrapping_add(prec_offset);
        for j in 0..iks_t() {
            let ks_ij = &ks_i[j];
            let k = (a_bar >> (32 - (j + 1) * basebit())) & ((1 << basebit()) - 1);
            if k != 0 {
                let ks_ijk = &ks_ij[k as usize];
                for x in 0..res.p.len() {
                    res.p[x] = res.p[x].wrapping_sub(ks_ijk.p[x]);
                }
            }
        }
    }

    return res;
}

pub fn generate_ksk(sKeyLv0: &Vec<u32>, sKeyLv1: &Vec<u32>) -> Vec<Vec<Vec<tlwe::TLWELv0>>> {
    let mut res: Vec<Vec<Vec<tlwe::TLWELv0>>> = Vec::new();

    for i in 0..N() {
        let mut ksk_i: Vec<Vec<tlwe::TLWELv0>> = Vec::new();
        for j in 0..iks_t() {
            let mut ksk_ij: Vec<tlwe::TLWELv0> = Vec::new();
            for k in 0..(1 << basebit()) {
                if k == 0 {
                    ksk_ij.push(tlwe::TLWELv0::new());
                    continue;
                }
                let p = ((k * sKeyLv1[i]) as f64) / ((1 << ((j + 1) * basebit())) as f64);
                ksk_ij.push(tlwe::tlweSymEncrypt(p, params::tlwe_lv0::ALPHA, sKeyLv0));
            }
            ksk_i.push(ksk_ij);
        }
        res.push(ksk_i);
    }

    return res;
}

pub fn hom_nand(
    tlwe_a: &tlwe::TLWELv0,
    tlwe_b: &tlwe::TLWELv0,
    ksk: &Vec<Vec<Vec<tlwe::TLWELv0>>>,
    b_key: &Vec<TRGSW>,
    test_vec: &trlwe::TRLWE,
    twist: &AlignedVec<c64>,
) -> tlwe::TLWELv0 {
    let mut tlwe_nand = -(tlwe_a + tlwe_b);
    *tlwe_nand.b_mut() = tlwe_nand.b() + utils::f64_to_torus(0.125);

    let trlwe = blind_rotate(&tlwe_nand, test_vec, b_key, twist);
    let tlwe_lv1 = trlwe::sample_extract_index(&trlwe, 0);
    let tlwe_bootstrapped = identity_key_switching(&tlwe_lv1, ksk);
    return tlwe_bootstrapped;
    //return tlwe_nand;
}

#[cfg(test)]
mod tests {
    use crate::mulfft;
    use crate::params;
    use crate::tlwe;
    use crate::trgsw::*;
    use crate::trlwe;
    use crate::utils;
    use rand::Rng;
    use std::time::{Duration, Instant};

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
            let h_u32 = utils::f64_to_torus_vec(&h);
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
        for i in 0..params::tlwe_lv0::N {
            keyLv0.push((rng.gen::<u8>() % 2) as u32);
        }
        for i in 0..N() {
            keyLv1.push((rng.gen::<u8>() % 2) as u32);
        }

        let mut bKey: Vec<TRGSW> = Vec::new();
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

            let tlwe = tlwe::tlweSymEncrypt(mu, params::tlwe_lv0::ALPHA, &keyLv0);
            let trlwe = blind_rotate(&tlwe, &test_vec, &bKey, &twist);
            let tlwe_lv1 = trlwe::sample_extract_index(&trlwe, 0);
            let dec = tlwe::tlweLv1SymDecrypt(&tlwe_lv1, &keyLv1);
            assert_eq!(plain_text, dec);
        }
    }

    #[test]
    fn test_identity_key_switching() {
        let mut rng = rand::thread_rng();
        let twist = mulfft::twist_gen(N());
        let mut keyLv0: Vec<u32> = Vec::new();
        let mut keyLv1: Vec<u32> = Vec::new();
        for i in 0..params::tlwe_lv0::N {
            keyLv0.push((rng.gen::<u8>() % 2) as u32);
        }
        for i in 0..N() {
            keyLv1.push((rng.gen::<u8>() % 2) as u32);
        }

        let ksk = generate_ksk(&keyLv0, &keyLv1);

        let try_num = 100;
        for i in 0..try_num {
            let plain_text = rng.gen::<u32>() % 2;
            let mut mu = 0.125;
            if plain_text == 0 {
                mu = -0.125;
            }

            let tlwe_lv1 = tlwe::tlweLv1SymEncrypt(mu, params::tlwe_lv0::ALPHA, &keyLv1);
            let tlwe_lv0 = identity_key_switching(&tlwe_lv1, &ksk);
            let dec = tlwe::tlweSymDecrypt(&tlwe_lv0, &keyLv0);
            assert_eq!(plain_text, dec);
        }
    }

    #[test]
    fn test_hom_nand() {
        let mut rng = rand::thread_rng();
        let twist = mulfft::twist_gen(N());
        let mut keyLv0: Vec<u32> = Vec::new();
        let mut keyLv1: Vec<u32> = Vec::new();
        for i in 0..params::tlwe_lv0::N {
            keyLv0.push((rng.gen::<u8>() % 2) as u32);
        }
        for i in 0..N() {
            keyLv1.push((rng.gen::<u8>() % 2) as u32);
        }

        let mut b_key: Vec<TRGSW> = Vec::new();
        for i in 0..keyLv0.len() {
            b_key.push(trgswSymEncrypt(keyLv0[i], trlwe::alpha(), &keyLv1, &twist));
        }

        let try_num = 100;
        let test_vec = generate_testvector();
        let ksk = generate_ksk(&keyLv0, &keyLv1);
        for i in 0..try_num {
            let plain_a = rng.gen::<u32>() % 2;
            let mut mu_a = 0.125;
            if plain_a == 0 {
                mu_a = -0.125;
            }
            let plain_b = rng.gen::<u32>() % 2;
            let mut mu_b = 0.125;
            if plain_b == 0 {
                mu_b = -0.125;
            }
            let mut nand = 0;
            if (plain_a & plain_b) == 0 {
                nand = 1;
            }

            let tlwe_a = tlwe::tlweSymEncrypt(mu_a, params::tlwe_lv0::ALPHA, &keyLv0);
            let tlwe_b = tlwe::tlweSymEncrypt(mu_b, params::tlwe_lv0::ALPHA, &keyLv0);
            let tlwe_nand = hom_nand(&tlwe_a, &tlwe_b, &ksk, &b_key, &test_vec, &twist);
            let dec = tlwe::tlweSymDecrypt(&tlwe_nand, &keyLv0);
            dbg!(plain_a);
            dbg!(plain_b);
            dbg!(nand);
            dbg!(dec);
            assert_eq!(nand, dec);
        }
    }

    #[test]
    fn test_hom_nand_bench() {
        let mut rng = rand::thread_rng();
        let twist = mulfft::twist_gen(N());
        let mut keyLv0: Vec<u32> = Vec::new();
        let mut keyLv1: Vec<u32> = Vec::new();
        for i in 0..params::tlwe_lv0::N {
            keyLv0.push((rng.gen::<u8>() % 2) as u32);
        }
        for i in 0..N() {
            keyLv1.push((rng.gen::<u8>() % 2) as u32);
        }

        let mut b_key: Vec<TRGSW> = Vec::new();
        for i in 0..keyLv0.len() {
            b_key.push(trgswSymEncrypt(keyLv0[i], trlwe::alpha(), &keyLv1, &twist));
        }

        let try_num = 10;
        let test_vec = generate_testvector();
        let ksk = generate_ksk(&keyLv0, &keyLv1);
        let plain_a = rng.gen::<u32>() % 2;
        let mut mu_a = 0.125;
        if plain_a == 0 {
            mu_a = -0.125;
        }
        let plain_b = rng.gen::<u32>() % 2;
        let mut mu_b = 0.125;
        if plain_b == 0 {
            mu_b = -0.125;
        }
        let mut nand = 0;
        if (plain_a & plain_b) == 0 {
            nand = 1;
        }

        let tlwe_a = tlwe::tlweSymEncrypt(mu_a, params::tlwe_lv0::ALPHA, &keyLv0);
        let tlwe_b = tlwe::tlweSymEncrypt(mu_b, params::tlwe_lv0::ALPHA, &keyLv0);
        let mut tlwe_nand = tlwe::TLWELv0::new();
        let start = Instant::now();
        for i in 0..try_num {
            tlwe_nand = hom_nand(&tlwe_a, &tlwe_b, &ksk, &b_key, &test_vec, &twist);
        }
        let end = start.elapsed();
        let exec_ms_per_gate = end.as_millis() as f64 / try_num as f64;
        println!("exec ms per gate : {} ms", exec_ms_per_gate);
        let dec = tlwe::tlweSymDecrypt(&tlwe_nand, &keyLv0);
        dbg!(plain_a);
        dbg!(plain_b);
        dbg!(nand);
        dbg!(dec);
        assert_eq!(nand, dec);
    }
}
