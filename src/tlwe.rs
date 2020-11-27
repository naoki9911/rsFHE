use crate::utils;
use crate::trlwe;
use rand::Rng;

pub struct TLWELv0 {
    pub p: Vec<u32>,
}

pub const fn n() -> i32 {
    630
}

pub const fn alpha() -> f64 {
    3.0517578125e-05
}

impl TLWELv0 {
    pub fn b(&self) -> u32 {
        return self.p[self.p.len()-1];
    }
}

pub fn tlweSymEncrypt(p: f64, alpha: f64, key: &Vec<u32>) -> TLWELv0 {
    let mut rng = rand::thread_rng();
    let mut tlwe: TLWELv0 = TLWELv0 { p: Vec::new() };
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

pub fn tlweSymDecrypt(tlwe: &TLWELv0, key: &Vec<u32>) -> u32 {
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

pub fn tlweLv1SymDecrypt(tlwe: &trlwe::TLWELv1, key :&Vec<u32>)  -> u32 {
    let mut tlwelv0 = TLWELv0{
        p :Vec::new()
    };
    for i in 0..tlwe.p.len() {
        tlwelv0.p.push(tlwe.p[i]);
    }
    return tlweSymDecrypt(&tlwelv0, key);
}

#[cfg(test)]
mod tests {
    use crate::tlwe::*;

    #[test]
    fn test_tlwe_enc_and_dec() {
        let mut rng = rand::thread_rng();

        // Generate 500bits secret key
        let mut key: Vec<u32> = Vec::new();
        let mut key_dirty: Vec<u32> = Vec::new();
        for i in 0..n() {
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
            let secret = tlweSymEncrypt(mu, alpha(), &key);
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
