use ark_ff::PrimeField;
use crate::multi_linear::MultiLinearPoly;

#[derive(Debug, Clone)]
pub struct Proof<F: PrimeField> {
    pub claimed_sums: Vec<F>,
    pub sum_polys: Vec<MultiLinearPoly<F>>,
    pub current_poly: MultiLinearPoly<F>,
}

pub struct ProverStruct<F: PrimeField> {
    // boolean hypercube computation
    // added the challenges array and final_univariate_poly to enable final testing
    pub bh_computation: MultiLinearPoly<F>,
    pub proof: Proof<F>,
}

impl<F: PrimeField> ProverStruct<F> {
    pub fn new(bh_computation: Vec<F>) -> Self {
        ProverStruct { 
            bh_computation: MultiLinearPoly::new(bh_computation.clone()), 
            proof: Proof { 
                claimed_sums: Vec::new(), 
                sum_polys: Vec::new(),
                current_poly: MultiLinearPoly::new(bh_computation)
            } 
        }
    }

    pub fn generate_proof(&mut self) -> Vec<(F, MultiLinearPoly<F>)> {
        let claimed_sum: F = self.proof.current_poly.computation.iter().sum();
        let half_len = self.proof.current_poly.computation.len() / 2;
    
        let (left, right) = self.proof.current_poly.computation.split_at(half_len);
        let left_sum: F = left.iter().sum();
        let right_sum = right.iter().sum();
    
        let sum_poly = MultiLinearPoly::new(vec![left_sum, right_sum]);
        self.proof.claimed_sums.push(claimed_sum);
        self.proof.sum_polys.push(sum_poly.clone());
    
        vec![(claimed_sum, sum_poly)]
    }

    pub fn next_poly(&mut self, challenge: F) {
        let new_poly = self.proof.current_poly.clone().partial_evaluate(challenge, 0);

        self.proof.current_poly = MultiLinearPoly::new(new_poly.computation.clone());
    }

    pub fn get_proof(&self) -> Proof<F> {
        self.proof.clone()
    }
}