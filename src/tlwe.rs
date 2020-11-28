use crate::params;
use crate::trlwe;
use crate::utils;
use rand::Rng;
use std::iter::Iterator;
use std::ops::Add;
use std::ops::Neg;

pub struct TLWELv0 {
    pub p: [u32; params::tlwe_lv0::N + 1],
}

impl TLWELv0 {
    pub fn new() -> TLWELv0 {
        return TLWELv0 {
            p: [0; params::tlwe_lv0::N + 1],
        };
    }

    pub fn b(&self) -> u32 {
        return self.p[params::tlwe_lv0::N];
    }

    pub fn b_mut(&mut self) -> &mut u32 {
        return &mut self.p[params::tlwe_lv0::N];
    }
}

impl Add for &TLWELv0 {
    type Output = TLWELv0;

    fn add(self, other: &TLWELv0) -> TLWELv0 {
        let mut res = TLWELv0::new();
        for ((rref, sval), oval) in res.p.iter_mut().zip(self.p.iter()).zip(other.p.iter()) {
            *rref = sval + oval;
        }
        return res;
    }
}

impl Neg for TLWELv0 {
    type Output = TLWELv0;

    fn neg(self) -> TLWELv0 {
        let mut res = TLWELv0::new();
        for (rref, sval) in res.p.iter_mut().zip(self.p.iter()) {
            *rref = 0u32.wrapping_sub(*sval);
        }

        return res;
    }
}

pub fn tlweSymEncrypt(p: f64, alpha: f64, key: &Vec<u32>) -> TLWELv0 {
    let mut rng = rand::thread_rng();
    let mut tlwe = TLWELv0::new();
    let mut inner_product: u32 = 0;

    for i in 0..key.len() {
        let rand_u32: u32 = rng.gen();
        inner_product = inner_product.wrapping_add(key[i] * rand_u32);
        tlwe.p[i] = rand_u32;
    }

    let mu = utils::f64_to_u32_torus(&vec![p]);
    let b = utils::gussian_32bit(&mu, alpha, 1);
    tlwe.p[params::tlwe_lv0::N] = inner_product.wrapping_add(b[0]);
    return tlwe;
}

pub fn tlweSymDecrypt(tlwe: &TLWELv0, key: &Vec<u32>) -> u32 {
    let mut inner_product: u32 = 0;
    for i in 0..key.len() {
        inner_product = inner_product.wrapping_add(tlwe.p[i] * key[i]);
    }

    let res_torus = (tlwe.p[params::tlwe_lv0::N].wrapping_sub(inner_product)) as i32;
    if res_torus < 0 {
        return 0;
    } else {
        return 1;
    }
}

pub fn tlweLv1SymDecrypt(tlwe: &trlwe::TLWELv1, key: &Vec<u32>) -> u32 {
    let mut inner_product: u32 = 0;
    for i in 0..key.len() {
        inner_product = inner_product.wrapping_add(tlwe.p[i] * key[i]);
    }

    let res_torus = (tlwe.p[key.len()].wrapping_sub(inner_product)) as i32;
    if res_torus < 0 {
        return 0;
    } else {
        return 1;
    }
}

pub fn tlweLv1SymEncrypt(p: f64, alpha: f64, key: &Vec<u32>) -> trlwe::TLWELv1 {
    let mut rng = rand::thread_rng();
    let mut tlwe: trlwe::TLWELv1 = trlwe::TLWELv1 { p: Vec::new() };
    let mut inner_product: u32 = 0;
    for i in 0..key.len() {
        let rand_u32: u32 = rng.gen();
        inner_product = inner_product.wrapping_add(key[i] * rand_u32);
        tlwe.p.push(rand_u32);
    }
    let mu = utils::f64_to_u32_torus(&vec![p]);
    let b = utils::gussian_32bit(&mu, alpha, 1);
    tlwe.p.push(inner_product.wrapping_add(b[0]));
    return tlwe;
}

#[cfg(test)]
mod tests {
    use crate::params;
    use crate::tlwe::*;

    #[test]
    fn test_tlwe_enc_and_dec() {
        let mut rng = rand::thread_rng();

        // Generate 500bits secret key
        let mut key: Vec<u32> = Vec::new();
        let mut key_dirty: Vec<u32> = Vec::new();
        for i in 0..params::tlwe_lv0::N {
            key.push((rng.gen::<u8>() % 2) as u32);
            key_dirty.push((rng.gen::<u8>() % 2) as u32);
        }

        let mut correct = 0;
        let try_num = 10000;

        for i in 0..try_num {
            let sample: u8 = rng.gen::<u8>() % 2;
            let mut mu = 0.125;
            if sample == 0 {
                mu = -0.125;
            }
            let secret = tlweSymEncrypt(mu, params::tlwe_lv0::ALPHA, &key);
            let plain = tlweSymDecrypt(&secret, &key) as u8;
            let plain_dirty = tlweSymDecrypt(&secret, &key_dirty) as u8;
            assert_eq!(plain, sample);
            if plain != plain_dirty {
                correct += 1;
            }
        }

        let probability = correct as f64 / try_num as f64;
        assert!(probability - 0.50 < 0.01);
    }
}
