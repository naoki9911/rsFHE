use fftw::array::AlignedVec;
use fftw::plan::*;
use fftw::types::*;
use rand::Rng;
use std::f64::consts::PI;

pub struct FFTPlan {
    pub forwardPlan: C2CPlan64,
    pub backwardPlan: C2CPlan64,
    pub twist: AlignedVec<c64>,
}

impl FFTPlan {
    pub fn new(N:usize) -> FFTPlan {
        return FFTPlan{
            forwardPlan: C2CPlan::aligned(&[N/2], Sign::Forward, Flag::MEASURE).unwrap(),
            backwardPlan: C2CPlan::aligned(&[N/2], Sign::Backward, Flag::MEASURE).unwrap(),
            twist:twist_gen(N),
        }
    }
}

pub fn twist_gen(N: usize) -> AlignedVec<c64> {
    let n: usize = N / 2;

    let mut array = AlignedVec::new(n);
    for k in 0..n {
        let value = (k as f64) * PI / (N as f64);
        // e^(i*value) = cos(value) + i*sin(value)
        array[k] = c64::new(value.cos(), value.sin());
    }

    return array;
}

pub fn twist_fft_1d(a: &Vec<f64>, plan: &mut FFTPlan) -> AlignedVec<c64> {
    let Ns = plan.twist.len();
    let mut t = AlignedVec::new(Ns);
    let mut res = AlignedVec::new(Ns);

    for i in 0..Ns {
        t[i] = plan.twist[i] * c64::new(a[i], a[Ns + i])
    }

    plan.forwardPlan.c2c(&mut t, &mut res).unwrap();

    return res;
}

pub fn twist_ifft_1d(input: &mut AlignedVec<c64>, plan: &mut FFTPlan) -> Vec<f64> {
    let Ns = plan.twist.len();
    let mut b = AlignedVec::new(Ns);
    let mut res: Vec<f64> = Vec::with_capacity(Ns * 2);

    plan.backwardPlan.c2c(input, &mut b).unwrap();
    for i in 0..Ns {
        b[i] = b[i] * plan.twist[i].conj();
        res.push(b[i].re / Ns as f64);
    }

    for i in 0..Ns {
        res.push(b[i].im / Ns as f64);
    }

    return res;
}

pub fn polynomial_mul(a: &Vec<i32>, b: &Vec<i32>, plan: &mut FFTPlan) -> Vec<u32> {
    let a_f64 = a.iter().map(|&e| e as f64).collect();
    let b_f64 = b.iter().map(|&e| e as f64).collect();
    let a_fft = twist_fft_1d(&a_f64, plan);
    let b_fft = twist_fft_1d(&b_f64, plan);

    let n = a_fft.len();
    let mut mul: AlignedVec<c64> = AlignedVec::new(n);
    for i in 0..n {
        mul[i] = a_fft[i] * b_fft[i];
    }

    let res = twist_ifft_1d(&mut mul, plan);
    return res
        .iter()
        .map(|&e| ((e.round() as i64) % 2i64.pow(32)) as u32)
        .collect();
}

pub fn polynomial_mul_u32(a: &Vec<u32>, b: &Vec<u32>, plan: &mut FFTPlan) -> Vec<u32> {
    let a_i32 = a.iter().map(|&e| e as i32).collect();
    let b_i32 = b.iter().map(|&e| e as i32).collect();
    return polynomial_mul(&a_i32, &b_i32, plan);
}

pub fn polynomial_mul_u32_1024(a: &[u32; 1024], b: &[u32; 1024], plan: &mut FFTPlan) -> [u32; 1024] {
    let a_f64 = a.iter().map(|&e| (e as i32) as f64).collect();
    let b_f64 = b.iter().map(|&e| (e as i32) as f64).collect();
    let a_fft = twist_fft_1d(&a_f64, plan);
    let b_fft = twist_fft_1d(&b_f64, plan);

    let n = a_fft.len();
    let mut mul: AlignedVec<c64> = AlignedVec::new(n);
    for i in 0..n {
        mul[i] = a_fft[i] * b_fft[i];
    }

    let res_f64 = twist_ifft_1d(&mut mul, plan);
    let mut res:[u32; 1024] = [0; 1024];
    let dividor:i64 = 2i64.pow(32);
    for (rref, f64val) in res.iter_mut().zip(res_f64.iter()) {
        *rref = ((f64val.round()) as i64 % dividor) as u32;
    }
    return res;
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

//pub fn polynomial_mul_lv2(a: &Vec<i64>, b: &Vec<i64>, twist: &AlignedVec<c64>) -> Vec<u64> {
//    let a_f64 = a.iter().map(|&e| e as f64).collect();
//    let b_f64 = b.iter().map(|&e| e as f64).collect();
//    let a_fft = twist_fft_1d(&a_f64, twist);
//    let b_fft = twist_fft_1d(&b_f64, twist);
//
//    let n = a_fft.len();
//    let mut mul:AlignedVec<c64> = AlignedVec::new(n);
//    for i in 0..n {
//        mul[i] = a_fft[i] * b_fft[i];
//    }
//
//    let res = twist_ifft_1d(&mut mul, twist);
//    return res.iter().map(|&e| ((e.round() as u128) % 2u128.pow(64)) as u64).collect();
//
//}

#[cfg(test)]
mod tests {
    use crate::mulfft::*;

    #[test]
    fn fft_add() {
        let a = [-3.0, -1.0, 1.0, 3.0];
        let b = [4.0, 3.0, 2.0, 1.0];
        let mut plan = FFTPlan::new(4);
        let a_fft = twist_fft_1d(&a.to_vec(), &mut plan);
        let b_fft = twist_fft_1d(&b.to_vec(), &mut plan);
        let mut sum: AlignedVec<c64> = AlignedVec::new(a_fft.len());
        for i in 0..a_fft.len() {
            sum[i] = a_fft[i] + b_fft[i];
        }
        let res = twist_ifft_1d(&mut sum, &mut plan);
        assert!(res[0] - 1.0 < 1e-10);
        assert!(res[1] - 2.0 < 1e-10);
        assert!(res[2] - 3.0 < 1e-10);
        assert!(res[3] - 4.0 < 1e-10);
    }

    #[test]
    fn fft_poly_mul() {
        let a: Vec<i32> = [-2, -1, 0, 1].to_vec();
        let b: Vec<i32> = [3, 4, 5, 6].to_vec();
        let mut plan = FFTPlan::new(4);
        let res = polynomial_mul(&a, &b, &mut plan);

        assert_eq!(res, vec![4294967292, 4294967280, 4294967276, 4294967282])
    }

    #[test]
    fn fft_ifft() {
        let mut a: Vec<u32> = Vec::new();
        let mut rng = rand::thread_rng();
        let mut plan = FFTPlan::new(1024);
        for i in 0..1024 {
            a.push(rng.gen::<u32>());
        }
        let a_f64 = a.iter().map(|&e| e as f64).collect();

        let mut a_fft = twist_fft_1d(&a_f64, &mut plan);
        let a_ifft = twist_ifft_1d(&mut a_fft, &mut plan);
        let res: Vec<u32> = a_ifft
            .iter()
            .map(|&e| ((e.round() as u64) % 2u64.pow(32)) as u32)
            .collect();
        for i in 0..1024 {
            assert_eq!(a[i], res[i]);
        }
    }

    //#[test]
    //fn fft_poly_mul_lv2(){
    //    let a:Vec<i64> = [-2, -1, 0, 1].to_vec();
    //    let b:Vec<i64> = [3, 4, 5, 6].to_vec();
    //    let twist = twist_gen(4);
    //    let res = polynomial_mul_lv2(&a, &b, &twist);

    //    println!("{:?}", res);
    //    assert_eq!(res, vec![4294967292, 4294967280, 4294967276, 4294967282])
    //}
}
