use crate::key;
use crate::mulfft;
use crate::params;
use crate::tlwe;
use crate::trlwe;
use crate::utils;
use std::convert::TryInto;

#[derive(Debug, Copy, Clone)]
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

#[derive(Debug, Copy, Clone)]
pub struct TRGSWLv1FFT {
    trlwe_fft: [trlwe::TRLWELv1FFT; params::trgsw_lv1::L * 2],
}

impl TRGSWLv1FFT {
    pub fn new(trgsw: &TRGSWLv1, plan: &mut mulfft::FFTPlan) -> TRGSWLv1FFT {
        return TRGSWLv1FFT {
            trlwe_fft: trgsw
                .trlwe
                .iter()
                .map(|t| trlwe::TRLWELv1FFT::new(t, plan))
                .collect::<Vec<trlwe::TRLWELv1FFT>>()
                .try_into()
                .unwrap(),
        };
    }

    pub fn new_dummy() -> TRGSWLv1FFT {
        return TRGSWLv1FFT {
            trlwe_fft: [trlwe::TRLWELv1FFT::new_dummy(); params::trgsw_lv1::L * 2],
        };
    }
}

pub fn trgswSymEncrypt(
    p: u32,
    alpha: f64,
    key: &key::SecretKeyLv1,
    plan: &mut mulfft::FFTPlan,
) -> TRGSWLv1 {
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

pub fn external_product_with_fft(
    trgsw_fft: &TRGSWLv1FFT,
    trlwe: &trlwe::TRLWELv1,
    cloud_key: &key::CloudKey,
    plan: &mut mulfft::FFTPlan,
) -> trlwe::TRLWELv1 {
    let dec = decomposition(&trlwe, cloud_key);

    let mut out_a_fft = [0.0f64; 1024];
    let mut out_b_fft = [0.0f64; 1024];

    const L: usize = params::trgsw_lv1::L;
    for i in 0..L * 2 {
        let dec_fft = plan.spqlios.ifft_1024(&dec[i]);
        fma_in_fd_1024(&mut out_a_fft, &dec_fft, &trgsw_fft.trlwe_fft[i].a);
        fma_in_fd_1024(&mut out_b_fft, &dec_fft, &trgsw_fft.trlwe_fft[i].b);
    }

    return trlwe::TRLWELv1 {
        a: plan.spqlios.fft_1024(&out_a_fft),
        b: plan.spqlios.fft_1024(&out_b_fft),
    };
}

fn fma_in_fd_1024(res: &mut [f64; 1024], a: &[f64; 1024], b: &[f64; 1024]) {
    for i in 0..512 {
        res[i] = a[i + 512] * b[i + 512] - res[i];
        res[i] = a[i] * b[i] - res[i];
        res[i + 512] += a[i] * b[i + 512] + a[i + 512] * b[i];
    }
}

pub fn decomposition(
    trlwe: &trlwe::TRLWELv1,
    cloud_key: &key::CloudKey,
) -> [[u32; params::trgsw_lv1::N]; params::trgsw_lv1::L * 2] {
    let mut res = [[0u32; params::trgsw_lv1::N]; params::trgsw_lv1::L * 2];

    let offset = cloud_key.decomposition_offset;
    const BGBIT: u32 = params::trgsw_lv1::BGBIT;
    const mask: u32 = (1 << params::trgsw_lv1::BGBIT) - 1;
    const half_bg: u32 = 1 << (params::trgsw_lv1::BGBIT - 1);

    for j in 0..params::trgsw_lv1::N {
        let tmp0 = trlwe.a[j].wrapping_add(offset);
        let tmp1 = trlwe.b[j].wrapping_add(offset);
        for i in 0..params::trgsw_lv1::L {
            res[i][j] = ((tmp0 >> (32 - ((i as u32) + 1) * BGBIT)) & mask).wrapping_sub(half_bg);
        }
        for i in 0..params::trgsw_lv1::L {
            res[i + params::trgsw_lv1::L][j] =
                ((tmp1 >> (32 - ((i as u32) + 1) * BGBIT)) & mask).wrapping_sub(half_bg);
        }
    }

    return res;
}

// if cond == 0 then in1 else in2
pub fn cmux(
    in1: &trlwe::TRLWELv1,
    in2: &trlwe::TRLWELv1,
    cond: &TRGSWLv1FFT,
    cloud_key: &key::CloudKey,
    plan: &mut mulfft::FFTPlan,
) -> trlwe::TRLWELv1 {
    let mut tmp = trlwe::TRLWELv1::new();
    const N: usize = params::trgsw_lv1::N;
    for i in 0..N {
        tmp.a[i] = in2.a[i].wrapping_sub(in1.a[i]);
        tmp.b[i] = in2.b[i].wrapping_sub(in1.b[i]);
    }

    let tmp2 = external_product_with_fft(cond, &tmp, cloud_key, plan);
    let mut res = trlwe::TRLWELv1::new();
    for i in 0..N {
        res.a[i] = tmp2.a[i].wrapping_add(in1.a[i]);
        res.b[i] = tmp2.b[i].wrapping_add(in1.b[i]);
    }

    return res;
}

pub fn blind_rotate(
    src: &tlwe::TLWELv0,
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
        res = cmux(
            &res,
            &res2,
            &cloud_key.bootstrapping_key[i as usize],
            cloud_key,
            plan,
        );
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
    cloud_key: &key::CloudKey,
    plan: &mut mulfft::FFTPlan,
) -> tlwe::TLWELv0 {
    let mut tlwe_nand = -(tlwe_a + tlwe_b);
    *tlwe_nand.b_mut() = tlwe_nand.b() + utils::f64_to_torus(0.125);

    let trlwe = blind_rotate(&tlwe_nand, cloud_key, plan);
    let tlwe_lv1 = trlwe::sample_extract_index(&trlwe, 0);
    let tlwe_bootstrapped = identity_key_switching(&tlwe_lv1, &cloud_key.key_switching_key);
    return tlwe_bootstrapped;
    //return tlwe_nand;
}

#[cfg(test)]
mod tests {
    use crate::key;
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
        let key = key::SecretKey::new();

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

            let c = trlwe::trlweSymEncrypt(
                &plain_text_enc,
                params::trlwe_lv1::ALPHA,
                &key.key_lv1,
                &mut plan,
            );
            let c_decomp = decomposition(&c, &cloud_key);
            let h_u32 = utils::f64_to_torus_vec(&h);
            let mut res = trlwe::TRLWELv1::new();
            for j in 0..N {
                let mut tmp0: u32 = 0;
                let mut tmp1: u32 = 0;
                for k in 0..params::trgsw_lv1::L {
                    tmp0 = tmp0.wrapping_add(c_decomp[k][j].wrapping_mul(h_u32[k]));
                    tmp1 = tmp1
                        .wrapping_add(c_decomp[k + params::trgsw_lv1::L][j].wrapping_mul(h_u32[k]));
                }
                res.a[j] = tmp0;
                res.b[j] = tmp1;
            }

            let dec = trlwe::trlweSymDecrypt(&res, &key.key_lv1, &mut plan);

            for j in 0..N {
                assert_eq!(plain_text[j], dec[j]);
            }
        }
    }

    #[test]
    fn test_external_product_with_fft() {
        const N: usize = params::trgsw_lv1::N;
        let mut rng = rand::thread_rng();
        let cloud_key = key::CloudKey::new_no_ksk();

        // Generate 1024bits secret key
        let key = key::SecretKey::new();

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

            let c = trlwe::trlweSymEncrypt(
                &plain_text_enc,
                params::trlwe_lv1::ALPHA,
                &key.key_lv1,
                &mut plan,
            );
            let p = trlwe::trlweSymDecrypt(&c, &key.key_lv1, &mut plan);
            let trgsw_true = trgswSymEncrypt(1, params::trgsw_lv1::ALPHA, &key.key_lv1, &mut plan);
            let trgsw_true_fft = TRGSWLv1FFT::new(&trgsw_true, &mut plan);
            let ext_c = external_product_with_fft(&trgsw_true_fft, &c, &cloud_key, &mut plan);
            let dec = trlwe::trlweSymDecrypt(&ext_c, &key.key_lv1, &mut plan);

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
        let key = key::SecretKey::new();
        let cloud_key = key::CloudKey::new_no_ksk();

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
            let c1 = trlwe::trlweSymEncrypt(&plain_text_enc_1, ALPHA, &key.key_lv1, &mut plan);
            let c2 = trlwe::trlweSymEncrypt(&plain_text_enc_2, ALPHA, &key.key_lv1, &mut plan);
            let trgsw_true = trgswSymEncrypt(1, ALPHA, &key.key_lv1, &mut plan);
            let trgsw_false = trgswSymEncrypt(0, ALPHA, &key.key_lv1, &mut plan);
            let trgsw_true_fft = TRGSWLv1FFT::new(&trgsw_true, &mut plan);
            let trgsw_false_fft = TRGSWLv1FFT::new(&trgsw_false, &mut plan);
            let enc_1 = cmux(&c1, &c2, &trgsw_false_fft, &cloud_key, &mut plan);
            let enc_2 = cmux(&c1, &c2, &trgsw_true_fft, &cloud_key, &mut plan);
            let dec_1 = trlwe::trlweSymDecrypt(&enc_1, &key.key_lv1, &mut plan);
            let dec_2 = trlwe::trlweSymDecrypt(&enc_2, &key.key_lv1, &mut plan);
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
        let key = key::SecretKey::new();
        let cloud_key = key::CloudKey::new(&key, &mut plan);

        let try_num = 10;
        for i in 0..try_num {
            let plain_text = rng.gen::<bool>();

            let tlwe =
                tlwe::TLWELv0::encrypt_bool(plain_text, params::tlwe_lv0::ALPHA, &key.key_lv0);
            let trlwe = blind_rotate(&tlwe, &cloud_key, &mut plan);
            let tlwe_lv1 = trlwe::sample_extract_index(&trlwe, 0);
            let dec = tlwe_lv1.decrypt_bool(&key.key_lv1);
            assert_eq!(plain_text, dec);
        }
    }

    #[test]
    fn test_identity_key_switching() {
        const N: usize = params::trgsw_lv1::N;
        let mut rng = rand::thread_rng();
        let key = key::SecretKey::new();
        let mut plan = mulfft::FFTPlan::new(N);
        let cloud_key = key::CloudKey::new(&key, &mut plan);

        let try_num = 100;
        for i in 0..try_num {
            let plain_text = rng.gen::<bool>();

            let tlwe_lv1 =
                tlwe::TLWELv1::encrypt_bool(plain_text, params::tlwe_lv1::ALPHA, &key.key_lv1);
            let tlwe_lv0 = identity_key_switching(&tlwe_lv1, &cloud_key.key_switching_key);
            let dec = tlwe_lv0.decrypt_bool(&key.key_lv0);
            assert_eq!(plain_text, dec);
        }
    }

    #[test]
    fn test_hom_nand() {
        const N: usize = params::trgsw_lv1::N;
        let mut rng = rand::thread_rng();
        let mut plan = mulfft::FFTPlan::new(N);
        let key = key::SecretKey::new();
        let cloud_key = key::CloudKey::new(&key, &mut plan);

        let try_num = 100;
        for i in 0..try_num {
            let plain_a = rng.gen::<bool>();
            let plain_b = rng.gen::<bool>();
            let nand = !(plain_a & plain_b);

            let tlwe_a =
                tlwe::TLWELv0::encrypt_bool(plain_a, params::tlwe_lv0::ALPHA, &key.key_lv0);
            let tlwe_b =
                tlwe::TLWELv0::encrypt_bool(plain_b, params::tlwe_lv0::ALPHA, &key.key_lv0);
            let tlwe_nand = hom_nand(&tlwe_a, &tlwe_b, &cloud_key, &mut plan);
            let dec = tlwe_nand.decrypt_bool(&key.key_lv0);
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
        let key = key::SecretKey::new();
        let cloud_key = key::CloudKey::new(&key, &mut plan);

        let mut b_key: Vec<TRGSWLv1> = Vec::new();
        for i in 0..key.key_lv0.len() {
            b_key.push(trgswSymEncrypt(
                key.key_lv0[i],
                params::trgsw_lv1::ALPHA,
                &key.key_lv1,
                &mut plan,
            ));
        }

        let try_num = 100;
        let plain_a = rng.gen::<bool>();
        let plain_b = rng.gen::<bool>();
        let nand = !(plain_a & plain_b);

        let tlwe_a = tlwe::TLWELv0::encrypt_bool(plain_a, params::tlwe_lv0::ALPHA, &key.key_lv0);
        let tlwe_b = tlwe::TLWELv0::encrypt_bool(plain_b, params::tlwe_lv0::ALPHA, &key.key_lv0);
        let mut tlwe_nand = tlwe::TLWELv0::new();
        println!("Started bechmark");
        let start = Instant::now();
        for i in 0..try_num {
            tlwe_nand = hom_nand(&tlwe_a, &tlwe_b, &cloud_key, &mut plan);
        }
        let end = start.elapsed();
        let exec_ms_per_gate = end.as_millis() as f64 / try_num as f64;
        println!("exec ms per gate : {} ms", exec_ms_per_gate);
        let dec = tlwe_nand.decrypt_bool(&key.key_lv0);
        dbg!(plain_a);
        dbg!(plain_b);
        dbg!(nand);
        dbg!(dec);
        assert_eq!(nand, dec);
    }
}
