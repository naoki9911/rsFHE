pub type Torus = u32;

pub mod tlwe_lv0 {

    pub const N: usize = 630;
    pub const ALPHA: f64 = 3.0517578125e-05;
}

pub mod tlwe_lv1 {

    pub const N: usize = 1024;
    pub const ALPHA: f64 = 2.98023223876953125e-08;
}

pub mod trlwe_lv1 {
    use crate::params;

    pub const N: usize = params::tlwe_lv1::N;
    pub const ALPHA: f64 = params::tlwe_lv1::ALPHA;
}

pub mod trgsw_lv1 {
    use crate::params;

    pub const N: usize = params::tlwe_lv1::N;
    pub const ALPHA: f64 = params::tlwe_lv1::ALPHA;
}
