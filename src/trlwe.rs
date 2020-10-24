use fftw::array::AlignedVec;
use fftw::types::*;
use rand::Rng;
use crate::mulfft;
use crate::utils;

pub fn trlweSymEncrypt(p:&Vec<f64>, alpha:f64, key:&Vec<u32>, twist: &AlignedVec<c64>) -> (Vec<u32>, Vec<u32>){
    let mut rng = rand::thread_rng();
    let mut a:Vec<u32> = Vec::new();
    for i in 0..key.len() {
        a.push(rng.gen());
    }
    let mut b = utils::gussian_32bit(&utils::f64_to_u32_torus(p), alpha, key.len());
    let a_i32 = a.iter().map(|&e| e as i32).collect();
    let key_i32 = key.iter().map(|&e| e as i32).collect();
    let poly_res = mulfft::polynomial_mul(&a_i32, &key_i32, twist);
    for i in 0..b.len() {
        b[i] = b[i].wrapping_add(poly_res[i]);
    }

    return (a, b);
}

pub fn trlweSymDecrypt(c:(&Vec<u32>, &Vec<u32>), key:&Vec<u32>, twist:&AlignedVec<c64>) -> Vec<u32> {
    let c_0_i32 = c.0.iter().map(|&e| e as i32).collect();
    let key_i32 = key.iter().map(|&e| e as i32).collect();
    let poly_res = mulfft::polynomial_mul(&c_0_i32, &key_i32, twist);
    let mut res:Vec<u32> = Vec::new();
    for i in 0..c.0.len() {
        let value = (c.1[i].wrapping_sub(poly_res[i])) as i32;
        if value < 0 {
            res.push(0);
        } else {
            res.push(1);
        }
    }
    return res;
}

#[cfg(test)]
mod tests {
    use crate::trlwe::*;
    use crate::mulfft;

    #[test]
    fn test_trlwe_enc_and_dec(){
        let mut rng = rand::thread_rng();

        // Generate 1024bits secret key
        let mut key:Vec<u32> = Vec::new();
        let mut key_dirty:Vec<u32> = Vec::new();
        for i in 0..1024 {
            key.push((rng.gen::<u8>() % 2) as u32);
            key_dirty.push((rng.gen::<u8>() % 2) as u32);
        }

        let alpha:f64 = 2.0f64.powf(-25.0);
        let twist = mulfft::twist_gen(1024);
        let mut correct = 0;
        let try_num = 500;

        for i in 0..try_num {
            let mut plain_text_enc:Vec<f64> = Vec::new();
            let mut plain_text:Vec<u32> = Vec::new();

            for j in 0..1024 {
                let sample:u32 = rng.gen::<u32>() % 2;
                let mut mu = 0.125;
                if sample == 0 {
                   mu = -0.125; 
                }
                plain_text.push(sample);
                plain_text_enc.push(mu);
            }

            let c = trlweSymEncrypt(&plain_text_enc, alpha, &key, &twist);
            let dec = trlweSymDecrypt((&c.0, &c.1), &key, &twist);
            let dec_dirty = trlweSymDecrypt((&c.0, &c.1), &key_dirty, &twist);

            for j in 0..1024 {
                assert_eq!(plain_text[j], dec[j]);
                if plain_text[j] != dec_dirty[j] {
                    correct += 1;
                }
            }
        }

        let probability = correct as f64 / (try_num * 1024) as f64;
        assert!(probability - 0.50 < 0.1);
    }
}