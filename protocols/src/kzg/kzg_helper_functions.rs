use crate::multi_linear::MultiLinearPoly;
use ark_ec::{pairing::Pairing, PrimeGroup};
use ark_ff::{PrimeField, Zero};

pub enum Operator {
    Add,
    Mul,
    Sub,
}

pub fn compute_commitment<F: PrimeField, P: Pairing>(
    poly: &MultiLinearPoly<F>,
    encrypted_basis: &[P::G1],
) -> P::G1 {
    let mut commitment = P::G1::zero();

    for (i, e_basis) in encrypted_basis.iter().enumerate() {
        commitment += e_basis.mul_bigint(poly.computation[i].into_bigint());
    }

    commitment
}

pub fn compute_poly_minus_v<F: PrimeField>(
    mut poly: MultiLinearPoly<F>,
    vars_to_open: &[F],
) -> MultiLinearPoly<F> {
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

pub fn compute_quotient<F: PrimeField>(poly: &MultiLinearPoly<F>) -> MultiLinearPoly<F> {
    let mid = poly.computation.len() / 2;

    let (first_half, second_half) = poly.computation.split_at(mid);
    let quotient = element_wise_op(second_half, first_half, Operator::Sub); // This is because Q(x) = arr of y_2 - arr of y_1

    MultiLinearPoly {
        computation: quotient,
    }
}

pub fn compute_remainder<F: PrimeField>(
    poly: MultiLinearPoly<F>,
    eval_point: F,
) -> MultiLinearPoly<F> {
    let mut this_poly = poly;
    let remainder = this_poly.partial_evaluate(eval_point, 0);

    remainder
}

pub fn blow_up<F: PrimeField>(
    poly: MultiLinearPoly<F>,
    blow_up_times: usize,
) -> MultiLinearPoly<F> {
    if !poly.computation.len().is_power_of_two() {
        panic!("The polynomial must have a power of 2 length");
    }
    let mut poly = poly;

    for _i in 0..blow_up_times {
        poly = MultiLinearPoly {
            computation: explode(poly),
        };
    }

    poly
}

fn explode<F: PrimeField>(poly: MultiLinearPoly<F>) -> Vec<F> {
    let n_bits = poly.computation.len();
    let total_combinations = 2 * n_bits;
    let mut blown_array = Vec::with_capacity(total_combinations);

    for _i in 0..2 {
        for val in poly.computation.iter() {
            blown_array.push(*val);
        }
    }

    blown_array
}

#[cfg(test)]
pub mod test {
    use crate::{kzg::trusted_setup::tests::setup, multi_linear::MultiLinearPoly};
    use ark_bls12_381::{Bls12_381, Fr as BlsFr, G1Affine};
    use ark_ec::AffineRepr;
    use ark_ff::PrimeField;

    pub fn poly_1() -> MultiLinearPoly<BlsFr> {
        // f(x) = 3ab + 4c
        MultiLinearPoly::new(&vec![
            BlsFr::from(0),
            BlsFr::from(4),
            BlsFr::from(0),
            BlsFr::from(4),
            BlsFr::from(0),
            BlsFr::from(4),
            BlsFr::from(3),
            BlsFr::from(7),
        ])
    }

    #[test]
    fn test_compute_quotient() {
        let poly = poly_1();
        let result = super::compute_quotient(&poly);
        let expected = MultiLinearPoly::new(&vec![
            BlsFr::from(0),
            BlsFr::from(0),
            BlsFr::from(3),
            BlsFr::from(3),
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_compute_commitment() {
        let trusted_setup = setup();
        let g1_generator = G1Affine::generator();
        let poly = poly_1();

        let result = super::compute_commitment::<BlsFr, Bls12_381>(&poly, &trusted_setup.g1_arr);

        let commitment = g1_generator.mul_bigint(BlsFr::from(42).into_bigint());

        assert_eq!(result, commitment);
    }

    #[test]
    fn test_test_compute_poly_minus_v() {
        let poly = poly_1();
        let vars_to_open = vec![BlsFr::from(6), BlsFr::from(4), BlsFr::from(0)];
        let result = super::compute_poly_minus_v(poly, &vars_to_open);
        let new_poly = vec![
            BlsFr::from(-72),
            BlsFr::from(-68),
            BlsFr::from(-72),
            BlsFr::from(-68),
            BlsFr::from(-72),
            BlsFr::from(-68),
            BlsFr::from(-69),
            BlsFr::from(-65),
        ];
        assert_eq!(new_poly, result.computation);
    }

    #[test]
    fn test_compute_remainder() {
        let poly = poly_1();
        let eval_point = BlsFr::from(6);

        let vars_to_open = vec![BlsFr::from(6), BlsFr::from(4), BlsFr::from(0)];
        let new_poly = super::compute_poly_minus_v(poly, &vars_to_open);

        let result = super::compute_remainder(new_poly, eval_point);
        let remainder = vec![
            BlsFr::from(-72),
            BlsFr::from(-68),
            BlsFr::from(-54),
            BlsFr::from(-50),
        ];
        assert_eq!(remainder, result.computation);
    }

    #[test]
    fn test_blow_up() {
        let poly = MultiLinearPoly::new(&vec![BlsFr::from(3), BlsFr::from(4)]);
        let result = super::blow_up(poly, 1);
        let new_poly = vec![
            BlsFr::from(3),
            BlsFr::from(4),
            BlsFr::from(3),
            BlsFr::from(4),
        ];
        assert_eq!(new_poly, result.computation);
    }
}
