mod key;
mod mulfft;
mod params;
mod spqlios;
mod tlwe;
mod trgsw;
mod trlwe;
mod utils;
use rand::Rng;

fn main() {
    let mut fft_plan = mulfft::FFTPlan::new(1024);
    let mut rng = rand::thread_rng();
    let secret_key = key::SecretKey::new();
    let cloud_key = key::CloudKey::new(&secret_key, &mut fft_plan);

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

    let tlwe_a = tlwe::tlweSymEncrypt(mu_a, params::tlwe_lv0::ALPHA, &secret_key.key_lv0);
    let tlwe_b = tlwe::tlweSymEncrypt(mu_b, params::tlwe_lv0::ALPHA, &secret_key.key_lv0);
    let tlwe_nand = trgsw::hom_nand(&tlwe_a, &tlwe_b, &cloud_key, &mut fft_plan);
    let dec = tlwe::tlweSymDecrypt(&tlwe_nand, &secret_key.key_lv0);
    println!("inA:{} inB:{} NandResult:{}", plain_a, plain_b, dec);
}
