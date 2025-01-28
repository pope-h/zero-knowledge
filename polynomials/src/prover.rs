use ark_ff::PrimeField;
use crate::{multi_linear::MultiLinearPoly, transcript::Transcript};
pub struct ProverStruct<F: PrimeField> {
    // boolean hypercube computation
    // added the challenges array and final_univariate_poly to enable final testing
    pub bh_computation: MultiLinearPoly<F>,
    pub challenges: Vec<F>,
    pub final_univariate_poly: MultiLinearPoly<F>,
}

impl<F: PrimeField> ProverStruct<F> {
    pub fn new(bh_computation: Vec<F>) -> Self {
        ProverStruct { bh_computation: MultiLinearPoly::new(bh_computation), challenges: Vec::new(), final_univariate_poly: MultiLinearPoly::new(Vec::new()) }
    }

    pub fn call_partial_evaluate(&mut self, eval_value: F, eval_value_position: usize) -> MultiLinearPoly<F> {
        MultiLinearPoly::partial_evaluate(&mut self.bh_computation, eval_value, eval_value_position)
    }

    pub fn convert_to_bytes(computation: Vec<F>) -> Vec<u8> {
        MultiLinearPoly::to_bytes(computation)
    }

    pub fn convert_from_bytes(bytes: &[u8]) -> Vec<F> {
        let converted = MultiLinearPoly::from_bytes(bytes);
        converted
    }

    pub fn generate_proof(&mut self) -> Vec<F> {
        let mut current_poly_ml = MultiLinearPoly::new(self.bh_computation.computation.clone());
        let mut transcript = Transcript::new();

        while current_poly_ml.computation.len() > 1 {
            println!("Starting computation is {:?}", current_poly_ml.computation);
            let claimed_sum: F = current_poly_ml.computation.iter().sum();

            let half_len = current_poly_ml.computation.len() / 2;
            let (left, right) = current_poly_ml.computation.split_at(half_len);
            let left_sum: F = left.iter().sum();
            let right_sum = right.iter().sum();

            let sum_poly = MultiLinearPoly::new(vec![left_sum, right_sum]);
            println!("Sum poly is {:?}", sum_poly);

            transcript.append(&ProverStruct::convert_to_bytes(vec![claimed_sum]));
            transcript.append(&ProverStruct::convert_to_bytes(sum_poly.computation.clone()));
            let challenge_bytes = transcript.challenge();
            let challenge = ProverStruct::convert_from_bytes(&challenge_bytes)[0];
            self.challenges.push(challenge);

            current_poly_ml = current_poly_ml.partial_evaluate(challenge, 0);
        }

        current_poly_ml.computation
    }

    pub fn verify_proof(&self) -> bool {
        let mut current_poly_ml = MultiLinearPoly::new(self.bh_computation.computation.clone());

        let final_output = current_poly_ml.evaluate(self.challenges.clone());
            
        final_output == self.final_univariate_poly
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    fn test_prover() {
        // 2a + 3b
        let computation = vec![Fq::from(0),
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
            Fq::from(12),];
        let mut prover = ProverStruct::new(computation);

        // let wrong_proof = vec![Fq::from(105)];
        // prover.final_univariate_poly = MultiLinearPoly::new(wrong_proof);

        let proof = prover.generate_proof();
        prover.final_univariate_poly = MultiLinearPoly::new(proof);
    
        assert!(prover.verify_proof(), "Proof verification failed");
    }
}