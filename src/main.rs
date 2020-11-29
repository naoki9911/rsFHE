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

    let plain_a = rng.gen::<bool>();
    let plain_b = rng.gen::<bool>();

    let tlwe_a = tlwe::TLWELv0::encrypt_bool(plain_a, params::tlwe_lv0::ALPHA, &secret_key.key_lv0);
    let tlwe_b = tlwe::TLWELv0::encrypt_bool(plain_b, params::tlwe_lv0::ALPHA, &secret_key.key_lv0);
    let tlwe_nand = trgsw::hom_nand(&tlwe_a, &tlwe_b, &cloud_key, &mut fft_plan);
    let dec = tlwe_nand.decrypt_bool(&secret_key.key_lv0);
    println!("inA:{} inB:{} NandResult:{}", plain_a, plain_b, dec);
}
