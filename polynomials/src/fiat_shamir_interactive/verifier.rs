// use ark_ff::PrimeField;
// use crate::{multi_linear::MultiLinearPoly, fiat_shamir_interactive::transcript::Transcript};

// pub struct VerifierStruct<F: PrimeField> {
//     pub bh_computation: MultiLinearPoly<F>,
//     pub challenges: Vec<F>,
//     pub final_eval_poly: Vec<F>,
// }

// impl<F: PrimeField> VerifierStruct<F> {
//     pub fn new(bh_computation: Vec<F>) -> Self {
//         VerifierStruct { bh_computation: MultiLinearPoly::new(bh_computation), challenges: Vec::new(), final_eval_poly: Vec::with_capacity(1) }
//     }

//     pub fn convert_to_bytes(computation: Vec<F>) -> Vec<u8> {
//         MultiLinearPoly::to_bytes(computation)
//     }

//     pub fn convert_from_bytes(bytes: &[u8]) -> Vec<F> {
//         let converted = MultiLinearPoly::from_bytes(bytes);
//         converted
//     }

//     pub fn variable_count(&self) -> u32 {
//         self.bh_computation.computation.len().ilog2()
//     }

//     pub fn check_proof(&mut self) {
//         let mut transcript = Transcript::new();
//         let count = self.variable_count();
//         let mut i = 0;

//         while i < (count - 1) {
//             let prover = ProverStruct::new(self.bh_computation.computation.clone());
//             let _ = prover.generate_proof();

//             let proof = prover.get_proof();
//             let claimed_sums = proof.claimed_sums;
//             let sum_polys = proof.sum_polys;

//             for i in 0..sum_polys.len() {
//                 let sum_poly_i = &sum_polys[i];
//                 if sum_poly_i.computation.iter().sum::<F>() != claimed_sums[i] {
//                     return;
//                 }

//                 transcript.append(&VerifierStruct::convert_to_bytes(vec![claimed_sums[i]]));
//                 transcript.append(&VerifierStruct::convert_to_bytes(sum_poly_i.computation.clone()));
//             }

//             let challenge_bytes = transcript.challenge();
//             let challenge = VerifierStruct::convert_from_bytes(&challenge_bytes)[0];
//             self.challenges.push(challenge);

//             let left = self.final_eval_poly.clone();
//             let right = self.final_eval_poly.clone();
//             let left_sum: F = left.iter().sum();
//             let right_sum = right.iter().sum();

//             let sum_poly = MultiLinearPoly::new(vec![left_sum, right_sum]);
//             self.final_eval_poly = sum_poly.computation.clone();
//             i += 1;
//         }
//     }

//     // pub fn check_proof(&mut self) -> bool {
//     //     let mut transcript = Transcript::new();

//     //     // call the get_proof function from the prover
//     //     let mut prover = ProverStruct::new(self.bh_computation.computation.clone());
//     //     let _ = prover.generate_proof();
//     //     let proof = prover.get_proof();

//     //     // get the claimed_sums and sum_polys from the proof
//     //     let claimed_sums = proof.claimed_sums;
//     //     let sum_polys = proof.sum_polys;

//     //     // add the sum_polys at an index to give the claimed_sum at that index
//     //     for i in 0..sum_polys.len() {
//     //         let sum_poly_i = &sum_polys[i];
//     //         if sum_poly_i.computation.iter().sum::<F>() != claimed_sums[i] {
//     //             return false;
//     //         }

//     //         transcript.append(&VerifierStruct::convert_to_bytes(vec![claimed_sums[i]]));
//     //         transcript.append(&VerifierStruct::convert_to_bytes(sum_poly_i.computation.clone()));
//     //         let challenge_bytes = transcript.challenge();
//     //         let challenge = VerifierStruct::convert_from_bytes(&challenge_bytes)[0];
//     //         self.challenges.push(challenge);

//     //         self.final_eval_poly = sum_poly_i.computation.clone();
//     //     }

//     //     true
//     // }
// }