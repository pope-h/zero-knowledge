use core::panic;

use crate::UnivariatePoly;
use rand::Rng;
use ark_ff::PrimeField;

#[derive(Debug)]
pub struct ShamirShare<F: PrimeField> {
    pub x: F,
    pub y: F,
}

impl<F: PrimeField> ShamirShare<F> {
    pub fn new(x: F, y: F) -> Self {
        ShamirShare { x, y }
    }
}

pub fn generate_shares<F: PrimeField>(secret: F, threshold: u8, num_shares: u8) -> Vec<ShamirShare<F>> {
    if threshold > num_shares {
        panic!("Threshold must be less than or equal to number of shares")
    }

    let polynomial = generate_random_polynomial(secret, threshold); // This is done filling up an array already having the secret at position 0 with other random numbers
    println!("Generated random polynomial: {:?}", polynomial);
    let x_values = generate_x_values(num_shares); // This is just generating numbers from 1 to be used as x values for the generated polynomial

    // This then pairs the x values with the y values gotten from evaluating the polynomial at the x values
    x_values
        .iter()
        .map(|x| ShamirShare::new(*x, polynomial.evaluate(*x)))
        .collect()
}

// generate x values i.e. 1, 2, 3, 4, 5, 6, 7
pub fn generate_x_values<F: PrimeField>(num_shares: u8) -> Vec<F> {
    (1..=num_shares).map(|x| F::from(x)).collect()
}

fn generate_random_polynomial<F: PrimeField>(secret: F, threshold: u8) -> UnivariatePoly<F> {
    let mut rng = rand::thread_rng();
    let mut coefficients = vec![secret];

    for _ in 1..threshold {
        coefficients.push(F::from(rng.gen_range(0, 100)));
    }

    UnivariatePoly::new(coefficients)
}

pub fn reconstruct_secret<F: PrimeField>(shares: &[ShamirShare<F>], threshold: u8) -> F {
    if shares.len() < threshold as usize {
        panic!("Not enough shares to reconstruct secret")
    }

    let xs: Vec<F> = shares.iter().map(|share| F::from(share.x)).collect();
    let ys: Vec<F> = shares.iter().map(|share| F::from(share.y)).collect();

    UnivariatePoly::interpolate(xs, ys).evaluate(F::from(0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    fn test_secret_sharing() {
        let secret = Fq::from(42);
        let threshold = 3;
        let num_shares = 7;

        // Generate shares
        let shares = generate_shares(secret, threshold, num_shares);
        println!("Generated {} shares:", shares.len());
        for share in &shares {
            println!("x: {}, y: {}", share.x, share.y);
        }

        // Test reconstruction with exactly threshold shares
        let reconstructed = reconstruct_secret(&shares[0..threshold as usize], threshold);
        assert_eq!(reconstructed, secret);

        // Test reconstruction with more than threshold shares
        let reconstructed = reconstruct_secret(&shares[0..4], threshold);
        assert_eq!(reconstructed, secret);
    }

    #[test]
    #[should_panic(expected = "Not enough shares")]
    fn test_insufficient_shares() {
        let shares = vec![ShamirShare::new(Fq::from(1), Fq::from(10)), ShamirShare::new(Fq::from(2), Fq::from(20))];
        reconstruct_secret(&shares, 3);
    }
}
