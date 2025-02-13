use crate::{multi_linear::MultiLinearPoly, transcript::Transcript};
use ark_ff::PrimeField;

#[derive(Debug, Clone)]
pub struct Proof<F: PrimeField> {
    pub init_poly_array: Vec<MultiLinearPoly<F>>,
    pub init_claimed_sum: F,
    pub sum_polys: Vec<MultiLinearPoly<F>>,
}

// The prover doesn't compute the claimed_sum in the proof fn but does it externally and passes it in to the proof fn
pub fn proof<F: PrimeField>(
    mut poly_array: Vec<MultiLinearPoly<F>>,
    init_claimed_sum: F,
) -> Proof<F> {
    if poly_array.len() != 2 {
        panic!("The number of polynomials must be 2");
    }

    if poly_array[0].computation.len() != poly_array[1].computation.len() {
        panic!("The polynomials must have the same length");
    }

    let init_poly_a = poly_array[0].computation.clone();
    let init_poly_b = poly_array[1].computation.clone();
    let mut transcript = Transcript::new();
    transcript.absorb(&MultiLinearPoly::to_bytes(
        poly_array[0].computation.clone(),
    ));
    transcript.absorb(&MultiLinearPoly::to_bytes(
        poly_array[1].computation.clone(),
    ));

    // let init_claimed_sum = poly.computation.iter().sum();
    let mut sum_polys = vec![];

    while poly_array[0].computation.len() > 1 {
        let half_len = poly_array[0].computation.len() / 2;

        let (left_a, right_a) = poly_array[0].computation.split_at(half_len);
        let (left_b, right_b) = poly_array[1].computation.split_at(half_len);

        let left_a_sum: F = left_a.iter().sum();
        let right_a_sum = right_a.iter().sum();
        let left_b_sum: F = left_b.iter().sum();
        let right_b_sum = right_b.iter().sum();

        let claimed_sum_a: F = poly_array[0].computation.iter().sum();
        let claimed_sum_b: F = poly_array[1].computation.iter().sum();

        let claimed_sum = claimed_sum_a + claimed_sum_b;

        let sum_poly_a = MultiLinearPoly::new(vec![left_a_sum, right_a_sum]);
        let sum_poly_b = MultiLinearPoly::new(vec![left_b_sum, right_b_sum]);

        let left_sum = sum_poly_a.computation[0] + sum_poly_b.computation[0];
        let right_sum = sum_poly_a.computation[1] + sum_poly_b.computation[1];

        let sum_poly = MultiLinearPoly::new(vec![left_sum, right_sum]);

        sum_polys.push(MultiLinearPoly {
            computation: sum_poly.computation.clone(),
        });

        transcript.absorb(&MultiLinearPoly::to_bytes(vec![claimed_sum]));
        transcript.absorb(&MultiLinearPoly::to_bytes(sum_poly.computation.clone()));
        let challenge_bytes = transcript.squeeze();
        let challenge = F::from_be_bytes_mod_order(&challenge_bytes);

        poly_array[0] = poly_array[0].partial_evaluate(challenge, 0);
        poly_array[1] = poly_array[1].partial_evaluate(challenge, 0);
    }

    Proof {
        init_poly_array: vec![
            MultiLinearPoly {
                computation: init_poly_a,
            },
            MultiLinearPoly {
                computation: init_poly_b,
            },
        ],
        init_claimed_sum,
        sum_polys,
    }
}

pub fn verify<F: PrimeField>(mut proof: Proof<F>) -> bool {
    let mut transcript = Transcript::new();
    transcript.absorb(&MultiLinearPoly::to_bytes(
        proof.init_poly_array[0].computation.clone(),
    ));
    transcript.absorb(&MultiLinearPoly::to_bytes(
        proof.init_poly_array[1].computation.clone(),
    ));

    let mut claimed_sum: F = proof.init_claimed_sum;
    let mut challenges: Vec<F> = vec![];

    for sum_poly in proof.sum_polys.iter() {
        let poly_sum: F = sum_poly.computation.iter().sum();
        if claimed_sum != poly_sum {
            return false;
        }

        transcript.absorb(&MultiLinearPoly::to_bytes(vec![claimed_sum]));
        transcript.absorb(&MultiLinearPoly::to_bytes(sum_poly.computation.clone()));
        let challenge_bytes = transcript.squeeze();
        let challenge = F::from_be_bytes_mod_order(&challenge_bytes);
        challenges.push(challenge);

        // verifier uses the (y_1 + (y_2 - y_1) * challenge) to evaluate the polynomial
        claimed_sum = sum_poly.computation[0]
            + ((sum_poly.computation[1] - sum_poly.computation[0]) * challenge);
    }

    let eval_a = proof.init_poly_array[0].evaluate(challenges.clone());
    let eval_b = proof.init_poly_array[1].evaluate(challenges);
    let final_eval = eval_a.computation[0] + eval_b.computation[0];

    final_eval == claimed_sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    fn test_proof() {
        let computation_a = vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)];
        let computation_b = vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)];
        let poly_a = MultiLinearPoly::new(computation_a);
        let poly_b = MultiLinearPoly::new(computation_b);

        let init_sum_a: Fq = poly_a.computation.iter().sum::<Fq>();
        let init_sum_b: Fq = poly_b.computation.iter().sum::<Fq>();

        let poly_array = vec![poly_a, poly_b];
        let init_claimed_sum = init_sum_a + init_sum_b;

        let proof = proof(poly_array, init_claimed_sum);
        dbg!(&proof);
        assert!(verify(proof));
    }
}
