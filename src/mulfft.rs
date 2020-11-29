use crate::spqlios;

pub struct FFTPlan {
    pub spqlios: *mut spqlios::SpqliosImpl,
    pub n: usize,
}

impl FFTPlan {
    pub fn new(n: usize) -> FFTPlan {
        return FFTPlan {
            spqlios: unsafe { spqlios::Spqlios_new(n as i32) },
            n: n,
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

#[cfg(test)]
mod tests {
    use crate::mulfft::*;
    use crate::params;
    use rand::Rng;

    #[test]
    fn test_spqlios_fft_ifft() {
        let n = 128;
        let mut plan = FFTPlan::new(n);
        let mut rng = rand::thread_rng();
        let mut a: Vec<u32> = vec![0u32; n];
        a.iter_mut().for_each(|e| *e = rng.gen::<u32>());

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
        let n = 128;
        let mut plan = FFTPlan::new(n);
        let mut rng = rand::thread_rng();
        let mut a: Vec<u32> = vec![0u32; n];
        let mut b: Vec<u32> = vec![0u32; n];
        a.iter_mut().for_each(|e| *e = rng.gen::<u32>());
        b.iter_mut()
            .for_each(|e| *e = rng.gen::<u32>() % params::trgsw_lv1::BG as u32);

        let spqlios_res = spqlios_poly_mul(&a, &b, &mut plan);
        let res = poly_mul(&a.to_vec(), &b.to_vec());
        for i in 0..n {
            let diff = res[i] as i32 - spqlios_res[i] as i32;
            assert!(diff < 2 && diff > -2);
        }
    }

    #[test]
    fn test_spqlios_fft_ifft_1024() {
        let mut plan = FFTPlan::new(1024);
        let mut rng = rand::thread_rng();
        let mut a = [0u32; 1024];
        a.iter_mut().for_each(|e| *e = rng.gen::<u32>());

        let a_fft = spqlios_ifft_1d_1024(&a, &mut plan);
        let res = spqlios_fft_1d_1024(&a_fft, &mut plan);
        for i in 0..1024 {
            let diff = a[i] as i32 - res[i] as i32;
            assert!(diff < 2 && diff > -2);
        }
    }

    #[test]
    fn test_spqlios_poly_mul_1024() {
        let mut plan = FFTPlan::new(1024);
        let mut rng = rand::thread_rng();
        let mut a = [0u32; 1024];
        let mut b = [0u32; 1024];
        a.iter_mut().for_each(|e| *e = rng.gen::<u32>());
        b.iter_mut()
            .for_each(|e| *e = rng.gen::<u32>() % params::trgsw_lv1::BG as u32);

        let spqlios_res = spqlios_poly_mul_1024(&a, &b, &mut plan);
        let res = poly_mul(&a.to_vec(), &b.to_vec());
        for i in 0..1024 {
            let diff = res[i] as i32 - spqlios_res[i] as i32;
            assert!(diff < 2 && diff > -2);
        }
    }

    fn poly_mul(a: &Vec<u32>, b: &Vec<u32>) -> Vec<u32> {
        let n = a.len();
        let mut res: Vec<u32> = vec![0u32; n];

        for i in 0..n {
            for j in 0..n {
                if i + j < n {
                    res[i + j] = res[i + j].wrapping_add(a[i].wrapping_mul(b[j]));
                } else {
                    res[i + j - n] = res[i + j - n].wrapping_sub(a[i].wrapping_mul(b[j]));
                }
            }
        }

        return res;
    }
    fn spqlios_ifft_1d(input: &Vec<u32>, plan: &mut FFTPlan) -> Vec<f64> {
        let mut res: Vec<f64> = vec![0.0f64; plan.n];

        unsafe {
            spqlios::Spqlios_ifft_lv1(plan.spqlios, res.as_mut_ptr(), input.as_ptr());
        }

        return res;
    }

    fn spqlios_fft_1d(input: &Vec<f64>, plan: &mut FFTPlan) -> Vec<u32> {
        let mut res: Vec<u32> = vec![0u32; plan.n];

        unsafe {
            spqlios::Spqlios_fft_lv1(plan.spqlios, res.as_mut_ptr(), input.as_ptr());
        }

        return res;
    }

    fn spqlios_poly_mul(a: &Vec<u32>, b: &Vec<u32>, plan: &mut FFTPlan) -> Vec<u32> {
        let a_ifft = spqlios_ifft_1d(a, plan);
        let b_ifft = spqlios_ifft_1d(b, plan);
        let mut mul = vec![0.0f64; plan.n];

        let ns = plan.n / 2;
        for i in 0..ns {
            let aimbim = a_ifft[i + ns] * b_ifft[i + ns];
            let arebim = a_ifft[i] * b_ifft[i + ns];
            mul[i] = a_ifft[i] * b_ifft[i] - aimbim;
            mul[i + ns] = a_ifft[i + ns] * b_ifft[i] + arebim;
        }

        return spqlios_fft_1d(&mul, plan);
    }
}
