use core::panic;

use rand::Rng;
use crate::UnivariatePoly;

#[derive(Debug)]
pub struct ShamirShare {
    pub x: f64,
    pub y: f64,
}

impl ShamirShare {
    pub fn new(x: f64, y: f64) -> Self {
        ShamirShare { x, y }
    }
    
}

pub fn generate_shares(secret: f64, threshold: u8, num_shares: u8) -> Vec<ShamirShare> {
    if threshold > num_shares {
        panic!("Threshold must be less than or equal to number of shares")
    }

    let polynomial = generate_random_polynomial(secret, threshold);
    println!("Generated random polynomial: {:?}", polynomial);
    let x_values = generate_x_values(num_shares);

    x_values.iter().map(|x| ShamirShare::new(*x, polynomial.evaluate(*x))).collect()
}

// generate random x values
pub fn generate_x_values(num_shares: u8) -> Vec<f64> {
    (1..=num_shares).map(|x| x as f64).collect()
}

fn generate_random_polynomial(secret: f64, threshold: u8) -> UnivariatePoly {
    let mut rng = rand::thread_rng();
    let mut coefficients = vec![secret];

    for _ in 1..threshold {
        coefficients.push(rng.gen_range(0, 100) as f64);
    }

    UnivariatePoly::new(coefficients)
}

pub fn reconstruct_secret(shares: &[ShamirShare], threshold: u8) -> f64 {
    if shares.len() < threshold as usize {
        panic!("Not enough shares to reconstruct secret")
    }

    let xs: Vec<f64> = shares.iter().map(|share| share.x).collect();
    let ys: Vec<f64> = shares.iter().map(|share| share.y).collect();

    UnivariatePoly::interpolate(xs, ys).evaluate(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_secret_sharing() {
        let secret = 42.0;
        let threshold = 3;
        let num_shares = 5;

        // Generate shares
        let shares = generate_shares(secret, threshold, num_shares);
        println!("Generated {} shares:", shares.len());
        for share in &shares {
            println!("x: {}, y: {}", share.x, share.y);
        }

        // Test reconstruction with exactly threshold shares
        let reconstructed = reconstruct_secret(&shares[0..threshold as usize], threshold);
        assert!((reconstructed - secret).abs() < 1e-10);

        // Test reconstruction with more than threshold shares
        let reconstructed = reconstruct_secret(&shares[0..4], threshold);
        assert!((reconstructed - secret).abs() < 1e-10);
    }

    #[test]
    #[should_panic(expected = "Not enough shares")]
    fn test_insufficient_shares() {
        let shares = vec![
            ShamirShare::new(1.0, 10.0),
            ShamirShare::new(2.0, 20.0),
        ];
        reconstruct_secret(&shares, 3);
    }
}