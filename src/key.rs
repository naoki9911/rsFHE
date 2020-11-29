use crate::params;
use crate::trlwe;
use crate::tlwe;
use crate::utils;

const TRGSWLV1_N:usize = params::trgsw_lv1::N;
const TRGSWLV1_IKS_T:usize = params::trgsw_lv1::IKS_T;
const TRGSWLV1_BASE:usize = 1 << params::trgsw_lv1::BASEBIT;
pub type KeySwitchingKey = Vec<tlwe::TLWELv0>;

pub struct CloudKey {
    pub decomposition_offset: u32,
    pub blind_rotate_testvec: trlwe::TRLWELv1,
    pub key_switching_key: KeySwitchingKey,
}

impl CloudKey {
    pub fn new(sKeyLv0:&Vec<u32>, sKeyLv1:&Vec<u32>) -> Self {
        return CloudKey {
            decomposition_offset: gen_decomposition_offset(),
            blind_rotate_testvec: gen_testvec(),
            key_switching_key: gen_key_switching_key(sKeyLv0, sKeyLv1),
        };
    }

    pub fn new_no_ksk() -> Self {
        return CloudKey {
            decomposition_offset: gen_decomposition_offset(),
            blind_rotate_testvec: gen_testvec(),
            key_switching_key: vec![tlwe::TLWELv0::new(); TRGSWLV1_BASE * TRGSWLV1_IKS_T * TRGSWLV1_N],
        };
    }
}

pub fn gen_decomposition_offset() -> u32 {
    let mut offset: u32 = 0;

    for i in 0..(params::trgsw_lv1::L as u32) {
        offset = offset.wrapping_add(
            params::trgsw_lv1::BG / 2 * (1 << (32 - (i + 1) * params::trgsw_lv1::BGBIT)),
        );
    }

    return offset;
}

pub fn gen_testvec() -> trlwe::TRLWELv1 {
    let mut testvec = trlwe::TRLWELv1::new();
    let b_torus = utils::f64_to_torus(0.125);
    for i in 0..params::trgsw_lv1::N {
        testvec.a[i] = 0;
        testvec.b[i] = b_torus;
    }

    return testvec;
}

pub fn gen_key_switching_key(sKeyLv0: &Vec<u32>, sKeyLv1: &Vec<u32>) -> KeySwitchingKey {
    const BASEBIT: usize = params::trgsw_lv1::BASEBIT;
    const IKS_T: usize = params::trgsw_lv1::IKS_T;
    let mut res  = vec![tlwe::TLWELv0::new(); TRGSWLV1_BASE * TRGSWLV1_IKS_T * TRGSWLV1_N];

    for i in 0..params::trgsw_lv1::N {
        for j in 0..IKS_T {
            for k in 0..TRGSWLV1_BASE {
                if k == 0 {
                    continue;
                }
                let p = ((k as u32 * sKeyLv1[i]) as f64) / ((1 << ((j + 1) * BASEBIT)) as f64);
                let idx = (TRGSWLV1_BASE * TRGSWLV1_IKS_T * i) + (TRGSWLV1_BASE * j) + k;
                res[idx] = tlwe::tlweSymEncrypt(p, params::tlwe_lv0::ALPHA, sKeyLv0);
            }
        }
    }

    return res;
}
