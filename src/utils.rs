use rand::thread_rng;
use rand_distr::{Normal, Distribution};

pub fn f64_to_u32_torus(d: &Vec<f64>) -> Vec<u32> {
    let mut res:Vec<u32> = Vec::new();
    for i in 0..d.len(){
        let torus = (d[i]%1.0) as f64 * 2u64.pow(32) as f64;
        res.push((torus as i64)as u32);
    }
    return res;
}

pub fn gussian_32bit(mu: &Vec<u32>, alpha: f64, size:usize) -> Vec<u32> {
    let normal = Normal::new(0.0, alpha).unwrap();
    let mut vec:Vec<u32> = Vec::new();
    for i in 0..size {
        let sample = normal.sample(&mut rand::thread_rng());
        vec.push(f64_to_u32_torus(&vec![sample])[0].wrapping_add(mu[i]));
    }

    return vec;
}


#[cfg(test)]
mod tests {
    use crate::utils::*;

    #[test]
    fn test_double_to_torust_32bit() {
        let torus = f64_to_u32_torus(&vec![3.141592]);
        assert_eq!(torus[0], 608133009);

        let torus2 = f64_to_u32_torus(&vec![2.71828]);
        assert_eq!(torus2[0], 3084989109);
    }

    #[test]
    fn test_gussian_32bit() {
        let torus = gussian_32bit(&vec![12], 0.1, 1);
        assert_eq!(torus.len(), 1);

        let torus2 = gussian_32bit(&vec![12,11], 0.1, 2);
        assert_eq!(torus2.len(), 2);
    }
}