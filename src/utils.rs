use rand::thread_rng;
use rand_distr::{Normal, Distribution};

pub fn double_to_torus_32bit(d: f64) -> u32 {
    let torus = (d%1.0) as f64 * 2u64.pow(32) as f64;
    return torus as u32;
}

pub fn gussian_32bit(mu: u32, alpha: f64, size:usize) -> Vec<u32> {
    let normal = Normal::new(0.0, alpha).unwrap();
    let mut vec:Vec<u32> = Vec::new();
    for i in 0..size {
        let sample = normal.sample(&mut rand::thread_rng());
        vec.push(double_to_torus_32bit(sample) + mu);
    }

    return vec;
}


#[cfg(test)]
mod tests {
    use crate::utils::*;

    #[test]
    fn test_double_to_torust_32bit() {
        let torus = double_to_torus_32bit(3.141592);
        assert_eq!(torus, 608133009);

        let torus2 = double_to_torus_32bit(2.71828);
        assert_eq!(torus2, 3084989109);
    }

    #[test]
    fn test_gussian_32bit() {
        let torus = gussian_32bit(12, 0.1, 1);
        assert_eq!(torus.len(), 1);

        let torus2 = gussian_32bit(12, 0.1, 2);
        assert_eq!(torus2.len(), 2);
    }
}