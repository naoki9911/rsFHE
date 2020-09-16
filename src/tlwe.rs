use rand::Rng;
use crate::utils;

pub fn tlweSymEncrypt(p:f64, alpha:f64, key:&Vec<u32>) -> Vec<u32> {
    let mut rng = rand::thread_rng();
    let mut a:Vec<u32> = Vec::new();
    let mut inner_product:u32 = 0;
    for i in 0..key.len() {
        let rand_u32: u32 = rng.gen();
        inner_product.wrapping_add(key[i] * rand_u32);
        a.push(rand_u32);
    }
    let mut b = utils::gussian_32bit(utils::f64_to_u32_torus(p), alpha, 1);
    a.push(inner_product + b[0]);
    return a;
}

pub fn tlweSymDecrypt(c:&Vec<u32>, key:&Vec<u32>) -> u32 {
    let mut inner_product:u32 = 0;
    for i in 0..key.len() {
        inner_product.wrapping_add(c[i] * key[i]);
    }

    let res_torus = (c[key.len()].wrapping_sub(inner_product)) as i32;
    if res_torus < 0 {
        return 0;
    } else {
        return 1;
    }
}

#[cfg(test)]
mod tests {
    use crate::tlwe::*;

    #[test]
    fn test_tlwe_enc_and_dec() {
        let key:Vec<u32> = vec![0, 1, 1, 0];
        let key_dirty:Vec<u32> = vec![1, 0, 0, 1];
        let alpha:f64 = 2.0f64.powf(-15.0);
        let mut rng = rand::thread_rng();
        for i in 0..100 {
            let sample:u8 = rng.gen::<u8>() % 2;
            let mut mu = 0.125;
            if sample == 0 {
               mu = -0.125; 
            }
            let secret = tlweSymEncrypt(mu, alpha, &key);
            let plain = tlweSymDecrypt(&secret, &key) as u8;
            let plain_dirty = tlweSymDecrypt(&secret, &key_dirty) as u8;
            println!("{} {}", plain, plain_dirty);
            assert_eq!(plain, sample);
            //assert_ne!(plain_dirty, sample)
        }
    }
}