use ark_ec::{pairing::Pairing, PrimeGroup};
use ark_ff::PrimeField;

#[derive(Debug)]
pub struct TrustedSetup<P: Pairing> {
    pub max_input: usize,
    pub g1_arr: Vec<P::G1>,
    pub g2_arr: Vec<P::G2>,
}

pub fn initialize<F: PrimeField, P: Pairing>(tau_arr: &[F]) -> TrustedSetup<P> {
    let max_arr_size = tau_arr.len();
    let lagrange_basis_arr = compute_lagrange_basis(&tau_arr);

    let g1_generator = P::G1::generator();
    let g2_generator = P::G2::generator();

    let encrypted_basis_poly = lagrange_basis_arr.iter().map(|val| g1_generator.mul_bigint(val.into_bigint())).collect();
    let encrypted_taus = tau_arr.iter().map(|tau| g2_generator.mul_bigint(tau.into_bigint())).collect();

    TrustedSetup { max_input: max_arr_size, g1_arr: encrypted_basis_poly, g2_arr: encrypted_taus }
}

// SHOULD IMPLEMENT THE CONTRIBUTE FUNCTION IN THE FUTURE
// SHOULD CHECK THAT THE INCOMING ARR IS SAME LEN AS THE MAX_INPUT

pub fn compute_lagrange_basis<F: PrimeField>(tau_arr: &[F]) -> Vec<F> {
    let poly_size = 2u32.pow(tau_arr.len() as u32) as usize;

    let num_bits = poly_size.trailing_zeros();

    let mut results = Vec::with_capacity(poly_size);

    for i in 0..poly_size {
        let mut product = F::one();

        // (0 & (1 << 5)) => 0 & (1 * 2^5) => 0 & 32 => 0
        for bit_position in 0..num_bits {
            // Calculate the bit position in Most-Significant-Bit-first order
            // normal order would go from LSB to MSB so we reverse it here
            // i.e. normal order: 000, 100, 010, 110, 001, 101, 011, 111
            // MSB-first order: 000, 001, 010, 011, 100, 101, 110, 111
            let msb_position = num_bits - 1 - bit_position;

            let bit_is_one = (i & (1 << msb_position)) != 0;    // Check if the bit at position bit_position is 1

            let val = if bit_is_one {
                tau_arr[bit_position as usize] // If bit is 1, use the variable directly
            } else {
                F::one() - tau_arr[bit_position as usize]  // If bit is 0, use (1 - variable)
            };

            product = product * val;
        }

        results.push(product);
    }

    results
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use ark_bls12_381::{Bls12_381, Fr as BlsFr, G1Affine};
    use ark_ec::AffineRepr;

    pub fn setup() -> TrustedSetup<Bls12_381> {
        let tau_arr = vec![BlsFr::from(5), BlsFr::from(2), BlsFr::from(3)];
        initialize::<BlsFr, Bls12_381>(&tau_arr)
    }

    #[test]
    fn test_compute_lagrange_basis() {
        let tau_arr = vec![BlsFr::from(5), BlsFr::from(2), BlsFr::from(3)];
        let result = compute_lagrange_basis(&tau_arr);
        let expected = vec![
            BlsFr::from(-8),
            BlsFr::from(12),
            BlsFr::from(16),
            BlsFr::from(-24),
            BlsFr::from(10),
            BlsFr::from(-15),
            BlsFr::from(-20),
            BlsFr::from(30),
        ];

        assert_eq!(result, expected);
    }

    #[test]
    fn test_initialize() {
        let tau_arr = vec![BlsFr::from(5), BlsFr::from(2), BlsFr::from(3)];
        let result = initialize::<BlsFr, Bls12_381>(&tau_arr);
        dbg!(&result);
    }

    #[test]
    fn test_negative() {
        let g1_generator = G1Affine::generator();
        let example = g1_generator.mul_bigint(BlsFr::from(-8).into_bigint());

        let a = example.mul_bigint(BlsFr::from(2).into_bigint());
        let b = g1_generator.mul_bigint(BlsFr::from(-16).into_bigint());
        // dbg!(&a, &b);   
        assert_eq!(a, b);    
    }
}