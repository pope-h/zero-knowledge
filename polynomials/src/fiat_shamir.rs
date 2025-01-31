use ark_ff::PrimeField;
use crate::{multi_linear::MultiLinearPoly, transcript::Transcript};

#[derive(Debug, Clone)]
pub struct Proof<F: PrimeField> {
    pub init_poly: MultiLinearPoly<F>,
    pub init_claimed_sum: F,
    pub sum_polys: Vec<MultiLinearPoly<F>>,
}


// The prover doesn't compute the claimed_sum in the proof fn but does it externally and passes it in to the proof fn
pub fn proof<F: PrimeField>(mut poly: MultiLinearPoly<F>, init_claimed_sum: F) -> Proof<F> {
    let init_poly = poly.computation.clone();
    let mut transcript = Transcript::new();
    transcript.absorb(&MultiLinearPoly::to_bytes(poly.computation.clone()));

    // let init_claimed_sum = poly.computation.iter().sum();
    let mut sum_polys = vec![];

    while poly.computation.len() > 1 {
        let half_len = poly.computation.len() / 2;
        let (left, right) = poly.computation.split_at(half_len);
        let left_sum: F = left.iter().sum();
        let right_sum = right.iter().sum();

        let claimed_sum: F = poly.computation.iter().sum();
        let sum_poly = MultiLinearPoly::new(vec![left_sum, right_sum]);
        // println!("Sum poly is {:?}", sum_poly);
        sum_polys.push(MultiLinearPoly { computation: sum_poly.computation.clone() });

        transcript.absorb(&MultiLinearPoly::to_bytes(vec![claimed_sum]));
        transcript.absorb(&MultiLinearPoly::to_bytes(sum_poly.computation.clone()));
        let challenge_bytes = transcript.squeeze();
        let challenge = MultiLinearPoly::from_bytes(&challenge_bytes)[0];

        poly = poly.partial_evaluate(challenge, 0);
    }

    Proof {
        init_poly: MultiLinearPoly { computation: init_poly },
        init_claimed_sum,
        sum_polys
    }
}

pub fn verify<F: PrimeField>(mut proof: Proof<F>) -> bool {
    let mut transcript = Transcript::new();

    // check that the polynomial is correct
    let verifier_claimed_sum = proof.init_poly.computation.iter().sum();
    if proof.init_claimed_sum != verifier_claimed_sum {
        return false;
    }
    transcript.absorb(&MultiLinearPoly::to_bytes(proof.init_poly.computation.clone()));

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
        let challenge = MultiLinearPoly::from_bytes(&challenge_bytes)[0];
        challenges.push(challenge);

        // verifier uses the (y_1 + (y_2 - y_1) * challenge) to evaluate the polynomial
        claimed_sum = sum_poly.computation[0] + ((sum_poly.computation[1] - sum_poly.computation[0]) * challenge);
    }

    let final_eval = proof.init_poly.evaluate(challenges);

    final_eval.computation[0] == claimed_sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fq;
    use crate::multi_linear::MultiLinearPoly;

    #[test]
    fn test_proof() {
        let poly = vec![MultiLinearPoly::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(3),
            Fq::from(5),
            Fq::from(9),
            Fq::from(8),
            Fq::from(12),
        ])];
        let init_claimed_sum = poly[0].computation.iter().sum();
        let proof = proof(poly[0].clone(), init_claimed_sum);
        let init_claimed_sum = proof.init_claimed_sum;

        println!("init_claimed_sum is {:?}", init_claimed_sum);
        assert_eq!(init_claimed_sum, Fq::from(48));
    }

    #[test]
    fn test_verify() {
        let poly = vec![MultiLinearPoly::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(3),
            Fq::from(5),
            Fq::from(9),
            Fq::from(8),
            Fq::from(12),
        ])];
        // let fake_claimed_sum = Fq::from(100);
        // let fake_proof = proof(poly[0].clone(), fake_claimed_sum);

        let init_claimed_sum = poly[0].computation.iter().sum();
        let proof = proof(poly[0].clone(), init_claimed_sum);
        let result = verify(proof);
        println!("Result is {:?}", result);
        assert!(result);
    }
}