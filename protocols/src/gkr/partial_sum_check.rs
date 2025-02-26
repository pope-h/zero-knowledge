use core::panic;

use crate::{
    gkr::product_poly::ProductPoly, multi_linear::MultiLinearPoly, transcript::Transcript,
    UnivariatePoly,
};
use ark_ff::PrimeField;

#[derive(Debug, Clone)]
pub struct Proof<F: PrimeField> {
    pub sum_poly: Vec<ProductPoly<F>>, // this is here so the verifier can compute the degree
    pub init_claimed_sum: F,
    pub challenges: Vec<F>,
    pub round_polys: Vec<Vec<F>>,
}

#[derive(Debug)]
pub struct SubClaim<F: PrimeField> {
    pub challenges: Vec<F>,
    pub last_claimed_sum: F,
}

// e.g. [[1, 2], [3, 4], [5, 6]] -> [[1 + 3 + 5], [2 + 4 + 6]] -> [9, 12]
pub fn reduce<F: PrimeField>(m_poly_array: Vec<Vec<F>>) -> Vec<F> {
    // let size = self.poly_array[0].computation.len() / 2;
    let mut new_array: Vec<F> = Vec::new();
    for i in 0..m_poly_array[0].len() {
        let sum = m_poly_array.iter().map(|array| array[i]).sum();

        new_array.push(sum);
    }

    new_array
}

pub fn proof<F: PrimeField>(mut sum_poly: Vec<ProductPoly<F>>, init_claimed_sum: F) -> Proof<F> {
    let mut initial_length = sum_poly[0].poly_array[0].computation.len().ilog2();
    let mut transcript = Transcript::new();
    let mut challenges: Vec<F> = vec![];

    let mut round_polys = vec![];

    while initial_length > 0 {
        let eval_array: Vec<Vec<F>> = sum_poly
            .iter()
            .map(|p_poly| p_poly.univariate_to_evaluation())
            .collect();
        let round_poly = reduce(eval_array);

        round_polys.push(round_poly.clone());
        transcript.absorb(&MultiLinearPoly::to_bytes(&round_poly));

        let challenge_bytes = transcript.squeeze();
        let challenge = F::from_be_bytes_mod_order(&challenge_bytes);
        challenges.push(challenge);

        sum_poly = sum_poly
            .iter()
            .map(|p_poly| p_poly.partial_evaluate(challenge, 0))
            .collect();

        initial_length -= 1;
    }

    Proof {
        sum_poly,
        init_claimed_sum,
        challenges,
        round_polys,
    }
}

// returns a struct of an array of challenges and last claimed_sum
pub fn verify<F: PrimeField>(proof: Proof<F>) -> SubClaim<F> {
    let mut transcript = Transcript::new();
    let mut claimed_sum: F = proof.init_claimed_sum;
    let mut challenges: Vec<F> = vec![];

    let degree = ProductPoly::get_degree(&proof.sum_poly[0]);
    let mut xs = Vec::with_capacity(degree);
    for i in 0..(degree + 1) {
        xs.push(F::from(i as u32));
    }

    for round_poly in proof.round_polys.iter() {
        let verifier_sum = round_poly[0] + round_poly[1];
        if claimed_sum != verifier_sum {
            panic!("Claimed sum does not match verifier sum");
        }

        transcript.absorb(&MultiLinearPoly::to_bytes(&round_poly));
        let challenge_bytes = transcript.squeeze();
        let challenge = F::from_be_bytes_mod_order(&challenge_bytes);
        challenges.push(challenge);

        let equation = UnivariatePoly::interpolate(&xs, &round_poly);
        claimed_sum = equation.evaluate(challenge);
    }

    SubClaim {
        challenges,
        last_claimed_sum: claimed_sum,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    fn test_proof() {
        // for a quick test, use [0, 0, 0, 2] and [0, 0, 0, 3]
        let poly_1 = MultiLinearPoly::new(vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ]);
        let poly_2 = MultiLinearPoly::new(vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ]);
        let prod_poly = ProductPoly::new(vec![poly_1, poly_2]);
        let init_claimed_sum = Fq::from(408);

        let proof = proof(vec![prod_poly.clone(), prod_poly], init_claimed_sum);
        dbg!(&proof);
        let verify = verify(proof);
        dbg!(&verify);
        assert_eq!(verify.challenges.len(), 3);
    }
}
