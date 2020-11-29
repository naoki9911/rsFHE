use crate::spqlios;
use rand::Rng;
use std::f64::consts::PI;

pub struct FFTPlan {
    pub spqlios: *mut spqlios::SpqliosImpl,
    pub n: usize,
}

impl FFTPlan {
    pub fn new(N: usize) -> FFTPlan {
        return FFTPlan {
            spqlios: unsafe { spqlios::Spqlios_new(N as i32) },
            n: N,
        };
    }
}

pub fn spqlios_ifft_1d_1024(input: &[u32; 1024], plan: &mut FFTPlan) -> [f64; 1024] {
    let src_const_ptr = input.as_ptr() as *const _;
    let mut res = Box::new([0.0f64; 1024]);
    let res_mut_ptr = Box::into_raw(res) as *mut _;

    unsafe {
        spqlios::Spqlios_ifft_lv1(plan.spqlios, res_mut_ptr, src_const_ptr);
    }
    res = unsafe { Box::from_raw(res_mut_ptr as *mut [f64; 1024]) };
    return *res;
}

pub fn spqlios_ifft_1d(input: &Vec<u32>, plan: &mut FFTPlan) -> Vec<f64> {
    let mut res: Vec<f64> = vec![0.0f64; plan.n];

    unsafe {
        spqlios::Spqlios_ifft_lv1(plan.spqlios, res.as_mut_ptr(), input.as_ptr());
    }

    return res;
}

pub fn spqlios_fft_1d_1024(input: &[f64; 1024], plan: &mut FFTPlan) -> [u32; 1024] {
    let src_const_ptr = input.as_ptr() as *const _;
    let mut res = Box::new([0u32; 1024]);
    let res_mut_ptr = Box::into_raw(res) as *mut _;

    unsafe {
        spqlios::Spqlios_fft_lv1(plan.spqlios, res_mut_ptr, src_const_ptr);
    }
    res = unsafe { Box::from_raw(res_mut_ptr as *mut [u32; 1024]) };
    return *res;
}

pub fn spqlios_fft_1d(input: &Vec<f64>, plan: &mut FFTPlan) -> Vec<u32> {
    let mut res: Vec<u32> = vec![0u32; plan.n];

    unsafe {
        spqlios::Spqlios_fft_lv1(plan.spqlios, res.as_mut_ptr(), input.as_ptr());
    }

    return res;
}

pub fn spqlios_poly_mul(a: &Vec<u32>, b: &Vec<u32>, plan: &mut FFTPlan) -> Vec<u32> {
    let a_ifft = spqlios_ifft_1d(a, plan);
    let b_ifft = spqlios_ifft_1d(b, plan);
    let mut mul = vec![0.0f64; plan.n];

    let Ns = plan.n / 2;
    for i in 0..Ns {
        let aimbim = a_ifft[i + Ns] * b_ifft[i + Ns];
        let arebim = a_ifft[i] * b_ifft[i + Ns];
        mul[i] = a_ifft[i] * b_ifft[i] - aimbim;
        mul[i + Ns] = a_ifft[i + Ns] * b_ifft[i] + arebim;
    }

    return spqlios_fft_1d(&mul, plan);
}

pub fn spqlios_poly_mul_1024(a: &[u32; 1024], b: &[u32; 1024], plan: &mut FFTPlan) -> [u32; 1024] {
    let a_ifft = spqlios_ifft_1d_1024(a, plan);
    let b_ifft = spqlios_ifft_1d_1024(b, plan);
    let mut mul = [0.0f64; 1024];

    for i in 0..512 {
        let aimbim = a_ifft[i + 512] * b_ifft[i + 512];
        let arebim = a_ifft[i] * b_ifft[i + 512];
        mul[i] = a_ifft[i] * b_ifft[i] - aimbim;
        mul[i + 512] = a_ifft[i + 512] * b_ifft[i] + arebim;
    }

    return spqlios_fft_1d_1024(&mul, plan);
}

pub fn poly_mul(a: &Vec<u32>, b: &Vec<u32>) -> Vec<u32> {
    let N = a.len();
    let mut res: Vec<u32> = Vec::new();

    for i in 0..N {
        res.push(0);
    }

    for i in 0..N {
        for j in 0..N {
            if (i + j < N) {
                res[i + j] = res[i + j].wrapping_add(a[i].wrapping_mul(b[j]));
            } else {
                res[i + j - N] = res[i + j - N].wrapping_sub(a[i].wrapping_mul(b[j]));
            }
        }
    }

    return res;
}

#[cfg(test)]
mod tests {
    use crate::mulfft::*;
    use crate::params;

    #[test]
    fn test_spqlios_fft_ifft() {
        let n = 128;
        let mut a: Vec<u32> = Vec::with_capacity(n);
        let mut rng = rand::thread_rng();
        let mut plan = FFTPlan::new(n);
        for i in 0..n {
            let tmp = rng.gen::<u32>();
            a.push(tmp);
        }

        let a_fft = spqlios_ifft_1d(&a, &mut plan);
        let res = spqlios_fft_1d(&a_fft, &mut plan);
        for i in 0..n {
            let diff = a[i] as i32 - res[i] as i32;
            assert!(diff < 2 && diff > -2);
            println!("{} {} {}", a_fft[i], a[i], res[i]);
        }
    }

    #[test]
    fn test_spqlios_poly_mul() {
        let n = 16;
        let mut a: Vec<u32> = Vec::with_capacity(n);
        let mut b: Vec<u32> = Vec::with_capacity(n);
        let mut rng = rand::thread_rng();
        let mut plan = FFTPlan::new(n);
        for i in 0..n {
            a.push(rng.gen::<u32>());
            b.push(rng.gen::<u32>() % params::trgsw_lv1::BG as u32);
        }

        let spqlios_res = spqlios_poly_mul(&a, &b, &mut plan);
        let res = poly_mul(&a.to_vec(), &b.to_vec());
        for i in 0..n {
            let diff = res[i] as i32 - spqlios_res[i] as i32;
            assert!(diff < 2 && diff > -2);
        }
    }

    #[test]
    fn test_spqlios_fft_ifft_1024() {
        let mut a = [0u32; 1024];
        let mut rng = rand::thread_rng();
        let mut plan = FFTPlan::new(1024);
        for i in 0..1024 {
            a[i] = rng.gen::<u32>();
        }

        let a_fft = spqlios_ifft_1d_1024(&a, &mut plan);
        let res = spqlios_fft_1d_1024(&a_fft, &mut plan);
        for i in 0..1024 {
            let diff = a[i] as i32 - res[i] as i32;
            assert!(diff < 2 && diff > -2);
        }
    }

    #[test]
    fn test_spqlios_poly_mul_1024() {
        let mut a = [0u32; 1024];
        let mut b = [0u32; 1024];
        let mut rng = rand::thread_rng();
        let mut plan = FFTPlan::new(1024);
        for i in 0..1024 {
            a[i] = rng.gen::<u32>();
            b[i] = rng.gen::<u32>() % params::trgsw_lv1::BG as u32;
        }

        let spqlios_res = spqlios_poly_mul_1024(&a, &b, &mut plan);
        let res = poly_mul(&a.to_vec(), &b.to_vec());
        for i in 0..1024 {
            let diff = res[i] as i32 - spqlios_res[i] as i32;
            assert!(diff < 2 && diff > -2);
        }
    }
}
