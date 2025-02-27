use ark_ec::{pairing::Pairing, PrimeGroup};
use ark_ff::{PrimeField, Zero};
use crate::multi_linear::MultiLinearPoly;

pub enum Operator {
    Add,
    Mul,
    Sub,
}

// vars_to_open ========== compute_binary_product
pub fn compute_lagrange_basis<F: PrimeField>(poly_size: usize, tau_arr: &[F]) -> Vec<F> {
    if !poly_size.is_power_of_two() {
        panic!("The polynomial must have a power of 2 number of coefficients");
    }

    let num_bits = poly_size.trailing_zeros();

    if tau_arr.len() != num_bits.try_into().unwrap() {
        panic!("The number of taus must be equal to the number of bits in the polynomial");
    }

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

pub fn compute_commitment<F: PrimeField, P: Pairing>(poly: MultiLinearPoly<F>, tau_arr: &[F]) -> P::G1 {
    let lagrange_basis: Vec<F> = compute_lagrange_basis(poly.computation.len(), tau_arr);
    let g1_generator = P::G1::generator();

    let encrypted_basis_poly: Vec<P::G1> = lagrange_basis.iter().map(|val| g1_generator.mul_bigint(val.into_bigint())).collect();
    let mut commitment = P::G1::zero();
    
    for (i, e_basis) in encrypted_basis_poly.iter().enumerate() {
        commitment += e_basis.mul_bigint(poly.computation[i].into_bigint());
    }

    let example = g1_generator.mul_bigint(F::from(-8).into_bigint());
    dbg!("First", example.mul_bigint(F::from(2).into_bigint()));
    dbg!("Second", g1_generator.mul_bigint(F::from(-16).into_bigint()));
    dbg!(g1_generator.mul_bigint(F::from(42).into_bigint()));
    commitment
}

pub fn compute_poly_sub_open_var<F: PrimeField>(mut poly: MultiLinearPoly<F>, vars_to_open: &[F]) -> MultiLinearPoly<F> {
    let eval_poly = poly.evaluate(vars_to_open);
    let v = eval_poly.computation[0];

    let sub_poly: Vec<F> = poly.computation.iter().map(|val| *val - v).collect();
    let result = MultiLinearPoly::new(&sub_poly);

    result
}

pub fn element_wise_op<F: PrimeField>(poly_a: &[F], poly_b: &[F], op: Operator) -> Vec<F> {
    if poly_a.len() != poly_b.len() {
        panic!("The two halves must have the same length");
    }

    let mut result = vec![F::zero(); poly_a.len()];

    for i in 0..poly_a.len() {
        result[i] = match op {
            Operator::Add => poly_a[i] + poly_b[i],
            Operator::Sub => poly_a[i] - poly_b[i],
            Operator::Mul => poly_a[i] * poly_b[i],
        };
    }

    result
}

pub fn compute_quotient<F: PrimeField>(poly: MultiLinearPoly<F>) -> Vec<F> {
    let mid = poly.computation.len() / 2;

    let (first_half, second_half) = poly.computation.split_at(mid);
    let quotient = element_wise_op(second_half, first_half, Operator::Sub); // This is because Q(x) = arr of y_2 - arr of y_1

    quotient
}

pub fn compute_remainder<F: PrimeField>(mut poly: MultiLinearPoly<F>, eval_point: F) -> MultiLinearPoly<F> {
    let remainder  = poly.partial_evaluate(eval_point, 0);

    remainder
}

#[cfg(test)]
mod test {
    use crate::{kzg_helper_functions::compute_lagrange_basis, multi_linear::MultiLinearPoly};
    use ark_bls12_381::Bls12_381;
    use ark_bn254::Fq;

    fn poly_1() -> MultiLinearPoly<Fq> {
        // f(x) = 3ab + 4c
        MultiLinearPoly::new(&vec![Fq::from(0), Fq::from(4), Fq::from(0), Fq::from(4), Fq::from(0), Fq::from(4), Fq::from(3), Fq::from(7)])
    }

    #[test]
    fn test_compute_binary_product() {
        let poly_size = poly_1().computation.len();
        let tau_arr = vec![Fq::from(5), Fq::from(2), Fq::from(3)];
        let result = compute_lagrange_basis(poly_size, &tau_arr);
        dbg!(&result);
    }

    #[test]
    fn test_compute_quotient() {
        let poly = poly_1();
        let result = super::compute_quotient(poly);
        dbg!(&result);
    }

    #[test]
    fn test_compute_commitment() {
        let poly = poly_1();
        let tau_arr = vec![Fq::from(5), Fq::from(2), Fq::from(3)];
        let result = super::compute_commitment::<Fq, Bls12_381>(poly, &tau_arr);
        dbg!(&result);
    }

    #[test]
    fn test_compute_remainder() {
        let poly = poly_1();
        let eval_point = Fq::from(6);
        let result = super::compute_remainder(poly, eval_point);
        dbg!(&result);
    }

    #[test]
    fn test_compute_poly_sub_open_var() {
        let poly = poly_1();
        let vars_to_open = vec![Fq::from(6), Fq::from(4), Fq::from(0)];
        let result = super::compute_poly_sub_open_var(poly, &vars_to_open);
        dbg!(Fq::from(-72));
        dbg!(&result);
    }
}