use ark_ff::PrimeField;
use crate::{multi_linear::MultiLinearPoly, prover::ProverStruct, transcript::Transcript};

pub struct VerifierStruct<F: PrimeField> {
    pub bh_computation: MultiLinearPoly<F>,
    pub challenges: Vec<F>,
    pub final_eval_poly: Vec<F>,
}

impl<F: PrimeField> VerifierStruct<F> {
    pub fn new(bh_computation: Vec<F>) -> Self {
        VerifierStruct { bh_computation: MultiLinearPoly::new(bh_computation), challenges: Vec::new(), final_eval_poly: Vec::with_capacity(1) }
    }

    pub fn convert_to_bytes(computation: Vec<F>) -> Vec<u8> {
        MultiLinearPoly::to_bytes(computation)
    }

    pub fn convert_from_bytes(bytes: &[u8]) -> Vec<F> {
        let converted = MultiLinearPoly::from_bytes(bytes);
        converted
    }

    pub fn check_proof(&mut self) -> bool {
        let mut transcript = Transcript::new();

        // call the get_proof function from the prover
        let mut prover = ProverStruct::new(self.bh_computation.computation.clone());
        let _ = prover.generate_proof();
        let proof = prover.get_proof();

        // get the claimed_sums and sum_polys from the proof
        let claimed_sums = proof.claimed_sums;
        let sum_polys = proof.sum_polys;

        // add the sum_polys at an index to give the claimed_sum at that index
        for i in 0..sum_polys.len() {
            let sum_poly_i = &sum_polys[i];
            if sum_poly_i.computation.iter().sum::<F>() != claimed_sums[i] {
                return false;
            }

            transcript.append(&VerifierStruct::convert_to_bytes(vec![claimed_sums[i]]));
            transcript.append(&VerifierStruct::convert_to_bytes(sum_poly_i.computation.clone()));
            let challenge_bytes = transcript.challenge();
            let challenge = VerifierStruct::convert_from_bytes(&challenge_bytes)[0];
            self.challenges.push(challenge);

            self.final_eval_poly = sum_poly_i.computation.clone();
        }

        true
    }

    // perform a final check to see if the final_eval_poly is equal to the final_eval_poly from the prover
    pub fn verify_proof(&self) -> bool {
        // convert the final_eval_poly to a univariate polynomial and evaluate it at the first challenge
        let mut final_eval: MultiLinearPoly<F> = MultiLinearPoly::new(self.final_eval_poly.clone());
        println!("Final eval is {:?}", final_eval);
        println!("Challenge is {:?}", self.challenges[self.challenges.len() - 1]);
        let final_eval_at_challenge = final_eval.partial_evaluate(self.challenges[self.challenges.len() - 1], 0);

        let mut this_computation = self.bh_computation.clone();
        let full_evaluation = MultiLinearPoly::new(vec![this_computation.evaluate(self.challenges.clone()).computation[0]]);
    
        final_eval_at_challenge == full_evaluation
    }
}

#[cfg(test)]
mod test {
    use crate::verifier::VerifierStruct;
    use ark_bn254::Fq;

    fn bh_computation() -> Vec<Fq> {
        vec![Fq::from(0),
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
            Fq::from(12)]
    }

    #[test]
    fn test_verify_proof() {
        let mut verifier = VerifierStruct::new(bh_computation());
        let check = verifier.check_proof();
        assert_eq!(check, true);

        let verify = verifier.verify_proof();
        assert_eq!(verify, true);
    }
}