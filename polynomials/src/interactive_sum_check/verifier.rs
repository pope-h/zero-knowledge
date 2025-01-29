use ark_ff::PrimeField;
use crate::{interactive_sum_check::transcript::Transcript, multi_linear::MultiLinearPoly};

pub struct VerifierStruct<F: PrimeField> {
    pub bh_computation: MultiLinearPoly<F>,
    pub challenges: Vec<F>,
    pub final_eval_poly: Vec<F>,
    pub transcript: Transcript,
}

impl<F: PrimeField> VerifierStruct<F> {
    pub fn new(bh_computation: Vec<F>) -> Self {
        VerifierStruct {
            bh_computation: MultiLinearPoly::new(bh_computation),
            challenges: Vec::new(),
            final_eval_poly: Vec::with_capacity(1),
            transcript: Transcript::new(),
        }
    }

    pub fn convert_to_bytes(computation: Vec<F>) -> Vec<u8> {
        MultiLinearPoly::to_bytes(computation)
    }

    pub fn convert_from_bytes(bytes: &[u8]) -> Vec<F> {
        let converted = MultiLinearPoly::from_bytes(bytes);
        converted
    }

    pub fn variable_count(&self) -> u32 {
        self.bh_computation.computation.len().ilog2()
    }

    pub fn check_proof(&mut self, proof_vec: Vec<(F, MultiLinearPoly<F>)>) -> bool {
        for (claimed_sum, sum_poly) in proof_vec {
            if sum_poly.computation.iter().sum::<F>() != claimed_sum {
                return false;
            }
            self.final_eval_poly = sum_poly.computation.clone();

            self.transcript.append(&VerifierStruct::convert_to_bytes(vec![claimed_sum]));
            self.transcript.append(&VerifierStruct::convert_to_bytes(sum_poly.computation.clone()));
        }

        true
    }

    pub fn generate_challenge(&mut self) -> F {
        let challenge_bytes = self.transcript.challenge();
        let challenge = VerifierStruct::convert_from_bytes(&challenge_bytes)[0];

        self.challenges.push(challenge);

        challenge
    }

    pub fn verify_proof(&mut self) -> bool {
        let mut final_eval: MultiLinearPoly<F> = MultiLinearPoly::new(self.final_eval_poly.clone());

        let final_eval_at_challenge = final_eval.partial_evaluate(self.challenges[self.challenges.len() - 1], 0);

        let mut this_computation = self.bh_computation.clone();
        let final_output = this_computation.evaluate(self.challenges.clone());

        final_output.computation[0] == final_eval_at_challenge.computation[0]
    }
}

/*
    STEPS FOR USING THE INTERACTIVE FIAT SHAMIR PROTOCOL
    1 => Create a new instance of the VerifierStruct
    2 => Call the generate_proof method on the ProverStruct instance
    3 => Call the check_proof method and if it returns true, call the generate_challenge method
    4 => Call the next_poly method on the ProverStruct instance
    5 => Repeat steps 2 to 4 for (variable_count.len() - 1) times
    6 => Repeat steps 2 & 3 (i.e. ensuring we don't send the last_challenge to the prover)
    7 => Finally call the verify_proof method on the VerifierStruct instance
 */

 #[cfg(test)]
mod test {
    use crate::interactive_sum_check::prover::ProverStruct;

    use super::*;
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
    let mut prover = ProverStruct::new(bh_computation());
    
    // Get the number of rounds needed (log2 of input size minus 1)
    let num_rounds = verifier.variable_count() - 1;
    
    for round in 0..num_rounds {
        let proof_array = prover.generate_proof();
        let verify = verifier.check_proof(proof_array);
        assert!(verify, "Proof verification failed at round {}", round);

        if verify {
            let challenge = verifier.generate_challenge();
            prover.next_poly(challenge);
        }
    }

    // Compute a final round that doesn't send challenge to the prover
    let final_proof = prover.generate_proof();
    if verifier.check_proof(final_proof) {
        verifier.generate_challenge();
        
        let final_eval = verifier.verify_proof();
        assert!(final_eval, "Final verification failed");

        // let final_eval_falsified = false;
        // assert!(final_eval_falsified, "Final verification failed");
    }

    assert!(verifier.verify_proof(), "Final verification failed");
    }
}