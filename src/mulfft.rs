use fftw::array::AlignedVec;
use fftw::types::*;
use fftw::plan::*;
use std::f64::consts::PI;

pub fn twist_gen(N: usize) -> AlignedVec<c64> {
    let n:usize = N / 2;

    let mut array = AlignedVec::new(n);
    for k in 0..n {
        let value = (k as f64)* PI / (N as f64);
        // e^(i*value) = cos(value) + i*sin(value)
        array[k] = c64::new(value.cos(), value.sin());
    }

    return array;
}


pub fn twist_fft_1d(a: &Vec<f64>, twist: &AlignedVec<c64>) -> AlignedVec<c64> {
    let Ns = twist.len();
    let mut t = AlignedVec::new(Ns);
    let mut res = AlignedVec::new(Ns);
    let mut plan: C2CPlan64 = C2CPlan::aligned(&[Ns], Sign::Forward, Flag::MEASURE).unwrap();

    for i in 0..Ns {
        t[i] = twist[i] * c64::new(a[i], a[Ns+i])
    }

    plan.c2c(&mut t, &mut res).unwrap();

    return  res;
}

pub fn twist_ifft_1d(input: &mut AlignedVec<c64>, twist: &AlignedVec<c64>) -> Vec<f64> {
    let Ns = twist.len();
    let mut b = AlignedVec::new(Ns);
    let mut res: Vec<f64> = Vec::with_capacity(Ns*2);
    let mut plan: C2CPlan64 = C2CPlan::aligned(&[Ns], Sign::Backward, Flag::MEASURE).unwrap();

    plan.c2c(input, &mut b).unwrap();
    for i in 0..Ns {
        b[i] = b[i] * twist[i].conj();
        res.push(b[i].re/Ns as f64);
    }

    for i in 0..Ns {
        res.push(b[i].im / Ns as f64);
    }

    return res;
}

#[cfg(test)]
mod tests {
    use crate::mulfft::*;

    #[test]
    fn fft_add(){
        let a = [-3.0,-1.0,1.0,3.0];
        let b = [4.0,3.0,2.0,1.0];
        let twist = twist_gen(4);
        let a_fft = twist_fft_1d(&a.to_vec(), &twist);
        let b_fft = twist_fft_1d(&b.to_vec(), &twist);
        let mut sum: AlignedVec<c64>= AlignedVec::new(a_fft.len());
        for i in 0..a_fft.len() {
            sum[i] = a_fft[i] + b_fft[i];
        }
        let res = twist_ifft_1d(&mut sum, &twist);
        println!("{:?}", res);
        assert!(res[0] - 1.0 < 1e-10);
        assert!(res[1] - 2.0 < 1e-10);
        assert!(res[2] - 3.0 < 1e-10);
        assert!(res[3] - 4.0 < 1e-10);
    }
}