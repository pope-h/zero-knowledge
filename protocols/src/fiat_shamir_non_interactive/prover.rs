use crate::{fiat_shamir_non_interactive::transcript::Transcript, multi_linear::MultiLinearPoly};
use ark_ff::PrimeField;

#[derive(Debug, Clone)]
pub struct Proof<F: PrimeField> {
    pub claimed_sums: Vec<F>,
    pub sum_polys: Vec<MultiLinearPoly<F>>,
}

pub struct FinalState<F: PrimeField> {
    pub challenges: Vec<F>,
    pub final_univariate_poly: Vec<F>,
}

pub struct ProverStruct<F: PrimeField> {
    // boolean hypercube computation
    // added the challenges array and final_univariate_poly to enable final testing
    pub bh_computation: MultiLinearPoly<F>,
    pub proof: Proof<F>,
    pub final_state: FinalState<F>,
}

impl<F: PrimeField> ProverStruct<F> {
    pub fn new(bh_computation: Vec<F>) -> Self {
        ProverStruct {
            bh_computation: MultiLinearPoly::new(bh_computation),
            proof: Proof {
                claimed_sums: Vec::new(),
                sum_polys: Vec::new(),
            },
            final_state: FinalState {
                challenges: Vec::new(),
                final_univariate_poly: Vec::new(),
            },
        }
    }

    pub fn convert_to_bytes(computation: Vec<F>) -> Vec<u8> {
        MultiLinearPoly::to_bytes(computation)
    }

    pub fn generate_proof(&mut self) -> Vec<F> {
        let mut current_poly_ml = MultiLinearPoly::new(self.bh_computation.computation.clone());
        let mut transcript = Transcript::new();
        transcript.append(&ProverStruct::convert_to_bytes(
            current_poly_ml.computation.clone(),
        ));

        while current_poly_ml.computation.len() > 1 {
            // println!("Starting computation is {:?}", current_poly_ml.computation);
            let claimed_sum: F = current_poly_ml.computation.iter().sum();

            let half_len = current_poly_ml.computation.len() / 2;
            let (left, right) = current_poly_ml.computation.split_at(half_len);
            let left_sum: F = left.iter().sum();
            let right_sum = right.iter().sum();

            let sum_poly = MultiLinearPoly::new(vec![left_sum, right_sum]);
            // println!("Sum poly is {:?}", sum_poly);
            self.proof.claimed_sums.push(claimed_sum);
            self.proof.sum_polys.push(MultiLinearPoly {
                computation: sum_poly.computation.clone(),
            });

            transcript.append(&ProverStruct::convert_to_bytes(vec![claimed_sum]));
            transcript.append(&ProverStruct::convert_to_bytes(
                sum_poly.computation.clone(),
            ));
            let challenge_bytes = transcript.challenge();
            let challenge = F::from_be_bytes_mod_order(&challenge_bytes);
            self.final_state.challenges.push(challenge);

            current_poly_ml = current_poly_ml.partial_evaluate(challenge, 0);
        }

        self.final_state.final_univariate_poly = current_poly_ml.computation.clone();
        current_poly_ml.computation
    }

    pub fn get_proof(&self) -> Proof<F> {
        // println!("Proof is {:?}", self.proof.clone());
        self.proof.clone()
    }

    pub fn verify_proof(&self) -> bool {
        let mut current_poly_ml = MultiLinearPoly::new(self.bh_computation.computation.clone());

        let final_output = current_poly_ml
            .evaluate(self.final_state.challenges.clone())
            .computation[0];

        final_output == self.final_state.final_univariate_poly[0]
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    fn test_prover() {
        // 2a + 3b
        let computation = vec![
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
        ];
        let mut prover = ProverStruct::new(computation);

        // let wrong_proof = vec![Fq::from(105)];
        // prover.final_state.final_univariate_poly = wrong_proof;

        let proof = prover.generate_proof();
        prover.final_state.final_univariate_poly = proof.clone();
        println!("PROOF is {:?}", proof);
        println!(
            "Final univariate poly is {:?}",
            prover.final_state.final_univariate_poly
        );

        assert!(prover.verify_proof(), "Proof verification failed");
    }
}
