// use ark_ff::PrimeField;
// use crate::{multi_linear::MultiLinearPoly, fiat_shamir_interactive::transcript::Transcript};

// #[derive(Debug, Clone)]
// pub struct Proof<F: PrimeField> {
//     pub claimed_sums: Vec<F>,
//     pub sum_polys: Vec<MultiLinearPoly<F>>,
// }

// pub struct ProverStruct<F: PrimeField> {
//     // boolean hypercube computation
//     // added the challenges array and final_univariate_poly to enable final testing
//     pub bh_computation: MultiLinearPoly<F>,
//     pub proof: Proof<F>,
// }

// impl<F: PrimeField> ProverStruct<F> {
//     pub fn new(bh_computation: Vec<F>) -> Self {
//         ProverStruct { bh_computation: MultiLinearPoly::new(bh_computation), proof: Proof { claimed_sums: Vec::new(), sum_polys: Vec::new() } }
//     }

//     pub fn call_partial_evaluate(&mut self, eval_value: F, eval_value_position: usize) -> MultiLinearPoly<F> {
//         MultiLinearPoly::partial_evaluate(&mut self.bh_computation, eval_value, eval_value_position)
//     }

//     pub fn generate_proof(&mut self) {
//         let mut current_poly = MultiLinearPoly::new(self.bh_computation.computation.clone());

//         while current_poly.computation.len() > 1 {
//             let claimed_sum: F = current_poly.computation.iter().sum();

//             let half_len = current_poly.computation.len() / 2;
//             let (left, right) = current_poly.computation.split_at(half_len);
//             let left_sum: F = left.iter().sum();
//             let right_sum = right.iter().sum();

//             let sum_poly = MultiLinearPoly::new(vec![left_sum, right_sum]);
//             self.proof.claimed_sums.push(claimed_sum);
//             self.proof.sum_polys.push(MultiLinearPoly { computation: sum_poly.computation.clone() });

//             get_challenge(&claimed_sum, &sum_poly);

//             current_poly = sum_poly;
//         }
//     }
// }