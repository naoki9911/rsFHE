use crate::key;
use crate::mulfft;
use crate::params;
use crate::tlwe;
use crate::trlwe;
use crate::utils;
use fftw::array::AlignedVec;
use fftw::types::*;
use std::convert::TryInto;

pub struct TRGSWLv1 {
    trlwe: [trlwe::TRLWELv1; params::trgsw_lv1::L * 2],
}

impl TRGSWLv1 {
    pub fn new() -> TRGSWLv1 {
        return TRGSWLv1 {
            trlwe: [trlwe::TRLWELv1::new(); params::trgsw_lv1::L * 2],
        };
    }
}

pub fn trgswSymEncrypt(p: u32, alpha: f64, key: &Vec<u32>, plan: &mut mulfft::FFTPlan) -> TRGSWLv1 {
    let mut p_f64: Vec<f64> = Vec::new();
    const L: usize = params::trgsw_lv1::L;
    for i in 0..L {
        let BgM: f64 = (params::trgsw_lv1::BG as f64).powf(((1 + i) as f64) * -1.0);
        p_f64.push(BgM);
    }
    let p_torus = utils::f64_to_torus_vec(&p_f64);

    let mut plain_zero: Vec<f64> = Vec::new();
    for i in 0..params::trgsw_lv1::N {
        plain_zero.push(0.0);
    }

    let mut trgsw = TRGSWLv1::new();
    trgsw
        .trlwe
        .iter_mut()
        .for_each(|e| *e = trlwe::trlweSymEncrypt(&plain_zero, alpha, key, plan));

    for i in 0..L {
        trgsw.trlwe[i].a[0] = trgsw.trlwe[i].a[0].wrapping_add(p * p_torus[i]);
        trgsw.trlwe[i + L].b[0] = trgsw.trlwe[i + L].b[0].wrapping_add(p * p_torus[i]);
    }
    return trgsw;
}

pub fn external_product(
    trgsw: &TRGSWLv1,
    trlwe: &trlwe::TRLWELv1,
    cloud_key: &key::CloudKey,
    plan: &mut mulfft::FFTPlan,
) -> trlwe::TRLWELv1 {
    let dec_a = decomposition(&trlwe.a, cloud_key);
    let dec_b = decomposition(&trlwe.b, cloud_key);
    let mut res = trlwe::TRLWELv1::new();

    const L: usize = params::trgsw_lv1::L;
    const N: usize = params::trgsw_lv1::N;
    for i in 0..L {
        let tmp = mulfft::polynomial_mul_u32_1024(&dec_a[i], &trgsw.trlwe[i].a, plan);
        for j in 0..N {
            res.a[j] = res.a[j].wrapping_add(tmp[j]);
        }

        let tmp = mulfft::spqlios_poly_mul_1024(&dec_b[i], &trgsw.trlwe[i + L].a, plan);
        for j in 0..N {
            res.a[j] = res.a[j].wrapping_add(tmp[j]);
        }

        let tmp = mulfft::spqlios_poly_mul_1024(&dec_a[i], &trgsw.trlwe[i].b, plan);
        for j in 0..N {
            res.b[j] = res.b[j].wrapping_add(tmp[j]);
        }

        let tmp = mulfft::spqlios_poly_mul_1024(&dec_b[i], &trgsw.trlwe[i + L].b, plan);
        for j in 0..N {
            res.b[j] = res.b[j].wrapping_add(tmp[j]);
        }
    }

    return res;
}

pub fn decomposition(
    a: &[u32; params::trgsw_lv1::N],
    cloud_key: &key::CloudKey,
) -> [[u32; params::trgsw_lv1::N]; params::trgsw_lv1::L] {
    let mut res = [[0u32; params::trgsw_lv1::N]; params::trgsw_lv1::L];
    let mut a_tilda = [0u32; params::trgsw_lv1::N];
    let offset = cloud_key.decomposition_offset;
    a_tilda
        .iter_mut()
        .zip(a.iter())
        .for_each(|(at_ref, a_val)| *at_ref = a_val.wrapping_add(offset));

    for i in 1..(params::trgsw_lv1::L + 1) as u32 {
        for j in 0..params::trgsw_lv1::N {
            let tmp = ((a_tilda[j as usize] >> (32 - params::trgsw_lv1::BGBIT * i))
                & (params::trgsw_lv1::BG - 1))
                .wrapping_sub(params::trgsw_lv1::BG / 2);
            res[(i - 1) as usize][j] = tmp;
        }
    }

    return res;
}

// if cond == 0 then in1 else in2
pub fn cmux(
    in1: &trlwe::TRLWELv1,
    in2: &trlwe::TRLWELv1,
    cond: &TRGSWLv1,
    cloud_key: &key::CloudKey,
    plan: &mut mulfft::FFTPlan,
) -> trlwe::TRLWELv1 {
    let mut tmp = trlwe::TRLWELv1::new();
    const N: usize = params::trgsw_lv1::N;
    for i in 0..N {
        tmp.a[i] = in2.a[i].wrapping_sub(in1.a[i]);
        tmp.b[i] = in2.b[i].wrapping_sub(in1.b[i]);
    }

    let tmp2 = external_product(cond, &tmp, cloud_key, plan);
    let mut res = trlwe::TRLWELv1::new();
    for i in 0..N {
        res.a[i] = tmp2.a[i].wrapping_add(in1.a[i]);
        res.b[i] = tmp2.b[i].wrapping_add(in1.b[i]);
    }

    return res;
}

pub fn blind_rotate(
    src: &tlwe::TLWELv0,
    bkey: &Vec<TRGSWLv1>,
    cloud_key: &key::CloudKey,
    plan: &mut mulfft::FFTPlan,
) -> trlwe::TRLWELv1 {
    const N: usize = params::trgsw_lv1::N;
    const NBIT: usize = params::trgsw_lv1::NBIT;
    let b_tilda = 2 * N - (((src.b() as usize) + (1 << (31 - NBIT - 1))) >> (32 - NBIT - 1));
    let mut res = trlwe::TRLWELv1 {
        a: poly_mul_with_X_k(&cloud_key.blind_rotate_testvec.a, b_tilda),
        b: poly_mul_with_X_k(&cloud_key.blind_rotate_testvec.b, b_tilda),
    };

    for i in 0..params::tlwe_lv0::N {
        let a_tilda =
            ((src.p[i as usize].wrapping_add((1 << (31 - NBIT - 1)))) >> (32 - NBIT - 1)) as usize;
        let res2 = trlwe::TRLWELv1 {
            a: poly_mul_with_X_k(&res.a, a_tilda),
            b: poly_mul_with_X_k(&res.b, a_tilda),
        };
        res = cmux(&res, &res2, &bkey[i as usize], cloud_key, plan);
    }

    return res;
}

pub fn poly_mul_with_X_k(a: &[u32; params::trgsw_lv1::N], k: usize) -> [u32; params::trgsw_lv1::N] {
    const N: usize = params::trgsw_lv1::N;

    let mut res: [u32; params::trgsw_lv1::N] = [0; params::trgsw_lv1::N];

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


pub fn identity_key_switching(
    src: &tlwe::TLWELv1,
    key_switching_key: &key::KeySwitchingKey,
) -> tlwe::TLWELv0 {
    const N: usize = params::trgsw_lv1::N;
    const BASEBIT: usize = params::trgsw_lv1::BASEBIT;
    const BASE: usize = 1 << BASEBIT;
    const IKS_T: usize = params::trgsw_lv1::IKS_T;
    let mut res = tlwe::TLWELv0::new();

    res.p[params::tlwe_lv0::N] = src.p[src.p.len() - 1];

    const prec_offset: u32 = 1 << (32 - (1 + BASEBIT * IKS_T));

    for i in 0..N {
        let a_bar = src.p[i].wrapping_add(prec_offset);
        for j in 0..IKS_T {
            let k = (a_bar >> (32 - (j + 1) * BASEBIT)) & ((1 << BASEBIT) - 1);
            if k != 0 {
                let idx = (BASE * IKS_T * i) + (BASE * j) + k as usize;
                for x in 0..res.p.len() {
                    res.p[x] = res.p[x].wrapping_sub(key_switching_key[idx].p[x]);
                }
            }
        }
    }

    return res;
}


pub fn hom_nand(
    tlwe_a: &tlwe::TLWELv0,
    tlwe_b: &tlwe::TLWELv0,
    b_key: &Vec<TRGSWLv1>,
    cloud_key: &key::CloudKey,
    plan: &mut mulfft::FFTPlan,
) -> tlwe::TLWELv0 {
    let mut tlwe_nand = -(tlwe_a + tlwe_b);
    *tlwe_nand.b_mut() = tlwe_nand.b() + utils::f64_to_torus(0.125);

    let trlwe = blind_rotate(&tlwe_nand, b_key, cloud_key, plan);
    let tlwe_lv1 = trlwe::sample_extract_index(&trlwe, 0);
    let tlwe_bootstrapped = identity_key_switching(&tlwe_lv1, &cloud_key.key_switching_key);
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
        const N: usize = params::trgsw_lv1::N;
        let mut rng = rand::thread_rng();
        let cloud_key = key::CloudKey::new_no_ksk();

        // Generate 1024bits secret key
        let mut key: Vec<u32> = Vec::new();
        for i in 0..N {
            key.push((rng.gen::<u8>() % 2) as u32);
        }

        let mut plan = mulfft::FFTPlan::new(N);
        let mut h: Vec<f64> = Vec::new();
        let try_num = 1000;

        for i in 1..params::trgsw_lv1::L + 1 {
            let tmp = (params::trgsw_lv1::BG as f64).powf(-(i as f64));
            h.push(tmp);
        }

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

            let c =
                trlwe::trlweSymEncrypt(&plain_text_enc, params::trlwe_lv1::ALPHA, &key, &mut plan);
            let c_decomp_1 = decomposition(&c.a, &cloud_key);
            let c_decomp_2 = decomposition(&c.b, &cloud_key);
            let h_u32 = utils::f64_to_torus_vec(&h);
            let mut res = trlwe::TRLWELv1::new();
            for j in 0..N {
                let mut tmp0: u32 = 0;
                let mut tmp1: u32 = 0;
                for k in 0..params::trgsw_lv1::L {
                    tmp0 = tmp0.wrapping_add(c_decomp_1[k][j].wrapping_mul(h_u32[k]));
                    tmp1 = tmp1.wrapping_add(c_decomp_2[k][j].wrapping_mul(h_u32[k]));
                }
                res.a[j] = tmp0;
                res.b[j] = tmp1;
            }

            let dec = trlwe::trlweSymDecrypt(&res, &key, &mut plan);

            for j in 0..N {
                assert_eq!(plain_text[j], dec[j]);
            }
        }
    }

    #[test]
    fn test_external_product() {
        const N: usize = params::trgsw_lv1::N;
        let mut rng = rand::thread_rng();
        let cloud_key = key::CloudKey::new_no_ksk();

        // Generate 1024bits secret key
        let mut key: Vec<u32> = Vec::new();
        for i in 0..N {
            key.push((rng.gen::<u8>() % 2) as u32);
        }

        let mut plan = mulfft::FFTPlan::new(1024);
        let try_num = 100;

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

            let c =
                trlwe::trlweSymEncrypt(&plain_text_enc, params::trlwe_lv1::ALPHA, &key, &mut plan);
            let p = trlwe::trlweSymDecrypt(&c, &key, &mut plan);
            let trgsw_true = trgswSymEncrypt(1, params::trgsw_lv1::ALPHA, &key, &mut plan);
            let ext_c = external_product(&trgsw_true, &c, &cloud_key, &mut plan);
            let dec = trlwe::trlweSymDecrypt(&ext_c, &key, &mut plan);

            for j in 0..N {
                assert_eq!(plain_text[j], p[j]);
            }
            for j in 0..N {
                assert_eq!(plain_text[j], dec[j]);
            }
        }
    }

    #[test]
    fn test_cmux() {
        const N: usize = params::trgsw_lv1::N;
        let mut rng = rand::thread_rng();
        let mut key: Vec<u32> = Vec::new();
        let cloud_key = key::CloudKey::new_no_ksk();
        for i in 0..N {
            key.push((rng.gen::<u8>() % 2) as u32);
        }

        let mut plan = mulfft::FFTPlan::new(N);
        let try_num = 100;
        for i in 0..try_num {
            let mut plain_text_enc_1: Vec<f64> = Vec::new();
            let mut plain_text_enc_2: Vec<f64> = Vec::new();
            let mut plain_text_1: Vec<u32> = Vec::new();
            let mut plain_text_2: Vec<u32> = Vec::new();

            for j in 0..N {
                let sample: u32 = rng.gen::<u32>() % 2;
                let mut mu = 0.125;
                if sample == 0 {
                    mu = -0.125;
                }
                plain_text_1.push(sample);
                plain_text_enc_1.push(mu);
            }
            for j in 0..N {
                let sample: u32 = rng.gen::<u32>() % 2;
                let mut mu = 0.125;
                if sample == 0 {
                    mu = -0.125;
                }
                plain_text_2.push(sample);
                plain_text_enc_2.push(mu);
            }
            const ALPHA: f64 = params::trgsw_lv1::ALPHA;
            let c1 = trlwe::trlweSymEncrypt(&plain_text_enc_1, ALPHA, &key, &mut plan);
            let c2 = trlwe::trlweSymEncrypt(&plain_text_enc_2, ALPHA, &key, &mut plan);
            let trgsw_true = trgswSymEncrypt(1, ALPHA, &key, &mut plan);
            let trgsw_false = trgswSymEncrypt(0, ALPHA, &key, &mut plan);
            let enc_1 = cmux(&c1, &c2, &trgsw_false, &cloud_key, &mut plan);
            let enc_2 = cmux(&c1, &c2, &trgsw_true, &cloud_key, &mut plan);
            let dec_1 = trlwe::trlweSymDecrypt(&enc_1, &key, &mut plan);
            let dec_2 = trlwe::trlweSymDecrypt(&enc_2, &key, &mut plan);
            for j in 0..N {
                assert_eq!(plain_text_1[j], dec_1[j]);
                assert_eq!(plain_text_2[j], dec_2[j]);
            }
        }
    }

    #[test]
    fn test_blind_rotate() {
        const N: usize = params::trgsw_lv1::N;
        let mut rng = rand::thread_rng();
        let mut plan = mulfft::FFTPlan::new(N);
        let mut keyLv0: Vec<u32> = Vec::new();
        let mut keyLv1: Vec<u32> = Vec::new();
        for i in 0..params::tlwe_lv0::N {
            keyLv0.push((rng.gen::<u8>() % 2) as u32);
        }
        for i in 0..N {
            keyLv1.push((rng.gen::<u8>() % 2) as u32);
        }
        let cloud_key = key::CloudKey::new(&keyLv0, &keyLv1);

        let mut bKey: Vec<TRGSWLv1> = Vec::new();
        for i in 0..keyLv0.len() {
            bKey.push(trgswSymEncrypt(
                keyLv0[i],
                params::trgsw_lv1::ALPHA,
                &keyLv1,
                &mut plan,
            ));
        }

        let try_num = 10;
        for i in 0..try_num {
            let plain_text = rng.gen::<u32>() % 2;
            let mut mu = 0.125;
            if plain_text == 0 {
                mu = -0.125;
            }

            let tlwe = tlwe::tlweSymEncrypt(mu, params::tlwe_lv0::ALPHA, &keyLv0);
            let trlwe = blind_rotate(&tlwe, &bKey, &cloud_key, &mut plan);
            let tlwe_lv1 = trlwe::sample_extract_index(&trlwe, 0);
            let dec = tlwe::tlweLv1SymDecrypt(&tlwe_lv1, &keyLv1);
            assert_eq!(plain_text, dec);
        }
    }

    #[test]
    fn test_identity_key_switching() {
        const N: usize = params::trgsw_lv1::N;
        let mut rng = rand::thread_rng();
        let twist = mulfft::twist_gen(N);
        let mut keyLv0: Vec<u32> = Vec::new();
        let mut keyLv1: Vec<u32> = Vec::new();
        for i in 0..params::tlwe_lv0::N {
            keyLv0.push((rng.gen::<u8>() % 2) as u32);
        }
        for i in 0..N {
            keyLv1.push((rng.gen::<u8>() % 2) as u32);
        }

        let cloud_key = key::CloudKey::new(&keyLv0, &keyLv1);

        let try_num = 100;
        for i in 0..try_num {
            let plain_text = rng.gen::<u32>() % 2;
            let mut mu = 0.125;
            if plain_text == 0 {
                mu = -0.125;
            }

            let tlwe_lv1 = tlwe::tlweLv1SymEncrypt(mu, params::tlwe_lv0::ALPHA, &keyLv1);
            let tlwe_lv0 = identity_key_switching(&tlwe_lv1, &cloud_key.key_switching_key);
            let dec = tlwe::tlweSymDecrypt(&tlwe_lv0, &keyLv0);
            assert_eq!(plain_text, dec);
        }
    }

    #[test]
    fn test_hom_nand() {
        const N: usize = params::trgsw_lv1::N;
        let mut rng = rand::thread_rng();
        let mut plan = mulfft::FFTPlan::new(N);
        let mut keyLv0: Vec<u32> = Vec::new();
        let mut keyLv1: Vec<u32> = Vec::new();
        for i in 0..params::tlwe_lv0::N {
            keyLv0.push((rng.gen::<u8>() % 2) as u32);
        }
        for i in 0..N {
            keyLv1.push((rng.gen::<u8>() % 2) as u32);
        }
        let cloud_key = key::CloudKey::new(&keyLv0, &keyLv1);

        let mut b_key: Vec<TRGSWLv1> = Vec::new();
        for i in 0..keyLv0.len() {
            b_key.push(trgswSymEncrypt(
                keyLv0[i],
                params::trgsw_lv1::ALPHA,
                &keyLv1,
                &mut plan,
            ));
        }

        let try_num = 100;
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
            let tlwe_nand = hom_nand(
                &tlwe_a, &tlwe_b, &b_key, &cloud_key, &mut plan,
            );
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
        const N: usize = params::trgsw_lv1::N;
        let mut rng = rand::thread_rng();
        let mut plan = mulfft::FFTPlan::new(N);
        let mut keyLv0: Vec<u32> = Vec::with_capacity(params::tlwe_lv0::N);
        let mut keyLv1: Vec<u32> = Vec::with_capacity(params::tlwe_lv1::N);
        for i in 0..params::tlwe_lv0::N {
            keyLv0.push((rng.gen::<u8>() % 2) as u32);
        }
        for i in 0..N {
            keyLv1.push((rng.gen::<u8>() % 2) as u32);
        }
        let cloud_key = key::CloudKey::new(&keyLv0, &keyLv1);

        let mut b_key: Vec<TRGSWLv1> = Vec::new();
        for i in 0..keyLv0.len() {
            b_key.push(trgswSymEncrypt(
                keyLv0[i],
                params::trgsw_lv1::ALPHA,
                &keyLv1,
                &mut plan,
            ));
        }

        let try_num = 100;
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
        println!("Started bechmark");
        let start = Instant::now();
        for i in 0..try_num {
            tlwe_nand = hom_nand(
                &tlwe_a, &tlwe_b, &b_key, &cloud_key, &mut plan,
            );
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
