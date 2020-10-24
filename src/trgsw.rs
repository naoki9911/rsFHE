//pub fn decomposition(trlwe:(&Vec<u32>, &Vec<u32>)) -> (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>) {
//    let offset = gen_offset();
//    println!("{}", offset);
//    let Bgbit:u32 = 10;
//    let halfBg:u32 = 1 << (Bgbit - 1);
//    let mask:u32 = (1 << Bgbit) - 1;
//    let N:u32 = 500;
//
//    let mut decomp0:Vec<u32> = Vec::new();
//    let mut decomp1:Vec<u32> = Vec::new();
//    let mut decomp2:Vec<u32> = Vec::new();
//    let mut decomp3:Vec<u32> = Vec::new();
//
//    for i in 0..N {
//        let temp0 = trlwe.0[i as usize].wrapping_add(offset);
//        let temp1 = trlwe.1[i as usize].wrapping_add(offset);
//
//        let mut tmp:u32 = ((temp0 >> (32 - (0 + 1) * Bgbit)) & mask).wrapping_sub(halfBg);
//        decomp0.push(tmp);
//        tmp = ((temp0 >> (32 - (1 + 1) * Bgbit)) & mask).wrapping_sub(halfBg);
//        decomp1.push(tmp);
//        tmp = ((temp1 >> (32 - (0 + 1) * Bgbit)) & mask).wrapping_sub(halfBg);
//        decomp2.push(tmp);
//        tmp = ((temp1 >> (32 - (1 + 1) * Bgbit)) & mask).wrapping_sub(halfBg);
//        decomp3.push(tmp);
//    }
//
//    return (decomp0, decomp1, decomp2, decomp3);
//}

pub fn decomposition(a:&Vec<u32>) -> (Vec<u32>, Vec<u32>) {
    let mut offset = 0;
    let Bgbit:u32 = 10;
    let Bg:u32 = 1 << Bgbit;
    let N:u32 = 1024;
    let offset = gen_offset();

    let mut a_tilda:Vec<u32> = Vec::new();
    let mut decomp0:Vec<u32> = Vec::new();
    let mut decomp1:Vec<u32> = Vec::new();

    for i in 0..N {
        a_tilda.push(a[i as usize].wrapping_add(offset));
    }

    for j in 0..N {
        let tmp0 = ((a_tilda[j as usize] >> (32-Bgbit*1))&(Bg-1)).wrapping_sub(Bg/2);
        decomp0.push(tmp0);
        let tmp1 = ((a_tilda[j as usize] >> (32-Bgbit*2))&(Bg-1)).wrapping_sub(Bg/2);
        decomp1.push(tmp1);
    }

    return (decomp0, decomp1);
}

pub fn gen_offset() -> u32 {
    let mut offset:u32 = 0;
    let Bgbit:u32 = 10;
    let Bg:u32 = 1 << Bgbit;

    for i in 0..2 {
        offset = offset.wrapping_add(Bg / 2 * (1 << (32 - (i + 1)*Bgbit)));
    }

    return offset;
}

#[cfg(test)]
mod tests {
    use crate::trlwe;
    use crate::utils;
    use crate::mulfft;
    use crate::trgsw::*;
    use rand::Rng;

    #[test]
    fn test_decomposition(){
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
        let mut h:Vec<f64> = Vec::new();
        let try_num = 1000;

        let Bgbit:u32 = 10;
        let Bg:u32 = 1 << Bgbit;
        for i in 0..2 {
            let mut tmp = (Bg as f64).powf(-(i+1) as f64);
            h.push(tmp);
        }

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


            let c = trlwe::trlweSymEncrypt(&plain_text_enc, alpha, &key, &twist);
            assert_eq!(c.0.len(), 1024);
            assert_eq!(c.1.len(), 1024);
            let c_decomp_1 = decomposition((&c.0));
            let c_decomp_2 = decomposition((&c.1));
            let h_u32 = utils::f64_to_u32_torus(&h);
            let mut rec0:Vec<u32> = Vec::new();
            let mut rec1:Vec<u32> = Vec::new();
            for j in 0..1024 {
                rec0.push((c_decomp_1.0[j].wrapping_mul(h_u32[0])).wrapping_add(c_decomp_1.1[j].wrapping_mul(h_u32[1])));
                rec1.push((c_decomp_2.0[j].wrapping_mul(h_u32[0])).wrapping_add(c_decomp_2.1[j].wrapping_mul(h_u32[1])));
            }

            let dec = trlwe::trlweSymDecrypt((&rec0, &rec1), &key, &twist);

            for j in 0..1024 {
                assert_eq!(plain_text[j], dec[j]);
                //if plain_text[j] != dec_dirty[j] {
                //}
            }
        }
    }
}