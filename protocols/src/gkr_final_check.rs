// use ark_ff::PrimeField;
// use crate::{gkr_protocol::{Circuit, GateOp}, multi_linear::MultiLinearPoly, partial_sum_check::{self, SubClaim}, product_poly::ProductPoly, transcript::Transcript};

// pub struct GKRProof<F: PrimeField> {
//     pub layer_proofs: Vec<LayerProof<F>>,
//     pub final_input_claims: Vec<F>
// }

// pub struct LayerProof<F: PrimeField> {
//     pub sum_check_rounds: Vec<Vec<F>>,  // polynomial coefficients for each round
//     pub final_w_claims: (F, F)  // (w_i+1(b*), w_i+1(c*)) claims for next layer
// }

// //GKRProof<F> 

// impl<F: PrimeField> Circuit<F> {
//     pub fn prove(&self, sub_claim: SubClaim<F>) {
//         let mut transcript = Transcript::new();
//         // let mut layer_proofs = Vec::new();

//         let r_b = sub_claim.challenges[0];
//         let r_c = sub_claim.challenges[1];
//         let last_claimed_sum = sub_claim.last_claimed_sum;
//         dbg!(r_b, r_c, last_claimed_sum);

//         // let circuit_inputs = self.inputs.clone();
//         let evaluated_circuit = self.evaluate();
//         dbg!(&evaluated_circuit);
//         let w_0 = evaluated_circuit[evaluated_circuit.len() -1].clone();
//         dbg!(&w_0);

//         let w_0_arr = if w_0.len() == 1 {
//             dbg!("single element");
//             vec![w_0[0], F::zero()]
//         } else if w_0.len().is_power_of_two() {
//             dbg!("power of two");
//             w_0
//         } else {
//             dbg!("not power of two");
//             let target_length = w_0.len().next_power_of_two();
//             let mut padded = w_0.clone();
//             padded.resize(target_length, F::zero());
//             padded
//         };
//         dbg!(&w_0_arr);

//         let m_0 = w_0_arr[0] + w_0_arr[1];  // claimed sum

//         transcript.absorb(&MultiLinearPoly::to_bytes(vec![m_0]));
//         let r_a = F::from_be_bytes_mod_order(&transcript.squeeze());

//         let mut w_0_mle = MultiLinearPoly::<F>::new(w_0_arr);
//         // m_0 = w_0(r)
//         // w_0(r) = f_0(r, b, c) => f_0(b, c)
//         let w_0_r = w_0_mle.partial_evaluate(r_a, 0);
//         dbg!(w_0_r);

//         // f_ri_b_c = [add_i_ri_b_c * (w_i+1_b + w_i+1_c)] + [mul_i_ri_b_c * (w_i+1_b * w_i+1_c)]
//         for layer_idx in (0..self.layers.len()).rev() {
//             // Get evaluations for current layer
//             let current_layer_eval = evaluated_circuit[layer_idx - 1].clone();
//             // let current_layer_eval = if layer_idx == self.layers.len() - 1 {
//             //     evaluated_circuit[layer_idx].clone()
//             // } else {
//             //     let prev_proof = &layer_proofs[layer_proofs.len() - 1];
//             //     vec![prev_proof.final_w_claims.0, prev_proof.final_w_claims.1]
//             // };

//             let n_bits = current_layer_eval.len().ilog2();

//             let (new_add, new_mul) = self.gkr_trick(r_b, r_c, layer_idx);
//             dbg!(&layer_idx);
//             dbg!(&new_add);
//             dbg!(&new_mul);

//             // [new_add(b, c) * (w(b) + w(c))] + [new_mul(b, c) * (w(b) * w(c))]
//             let (w_i_b_exploded, w_i_c_exploded) = Circuit::<F>::eval_layer_i(current_layer_eval, n_bits);
//             dbg!(&w_i_b_exploded);
//             dbg!(&w_i_c_exploded);
//             dbg!("got here");
//             let mut w_i_b_mle = MultiLinearPoly::<F>::new(w_i_b_exploded);
//             dbg!("check point");
//             let mut w_i_c_mle = MultiLinearPoly::<F>::new(w_i_c_exploded);

//             let w_i_b = w_i_b_mle.partial_evaluate(r_b, 0);
//             dbg!("checking here");
//             let w_i_c = w_i_c_mle.partial_evaluate(r_c, 0);

//             // f(b, c) using new_add and new_mul
//             // [new_add * (w_i_b + w_i_c)] + [new_mul * (w_i_b * w_i_c)]
//             let w_i_b_c_add = Circuit::<F>::element_wise_op(w_i_b.computation.clone(), w_i_c.computation.clone(), GateOp::Add);
//             let w_i_b_c_mul = Circuit::<F>::element_wise_op(w_i_b.computation, w_i_c.computation, GateOp::Mul);

//             let p_poly_1 = ProductPoly::new(vec![new_add, MultiLinearPoly { computation: w_i_b_c_add }]);
//             let p_poly_2 = ProductPoly::new(vec![new_mul, MultiLinearPoly { computation: w_i_b_c_mul }]);

//             let result = partial_sum_check::proof::<F>(vec![p_poly_1, p_poly_2], m_0);
//             dbg!(result);
//         }
//     }
// }

// #[cfg(test)]
// mod test {
//     use crate::gkr_protocol::{Gate, Layer};

//     use super::*;
//     use ark_bn254::Fq;

//     #[test]
//     fn test_gkr_proof() {
//         let inputs = vec![
//             Fq::from(1),
//             Fq::from(2),
//             Fq::from(3),
//             Fq::from(4),
//             Fq::from(5),
//             Fq::from(6),
//             Fq::from(7),
//             Fq::from(8),
//         ];
//         let mut circuit = Circuit::new(inputs);

//         let layer_1 = Layer {
//             gates: vec![
//                 Gate {
//                     left: 0,
//                     right: 1,
//                     op: GateOp::Add,
//                     output: 0,
//                 },
//                 Gate {
//                     left: 2,
//                     right: 3,
//                     op: GateOp::Mul,
//                     output: 1,
//                 },
//                 Gate {
//                     left: 4,
//                     right: 5,
//                     op: GateOp::Mul,
//                     output: 2,
//                 },
//                 Gate {
//                     left: 6,
//                     right: 7,
//                     op: GateOp::Mul,
//                     output: 3,
//                 },
//             ],
//         };

//         let layer_2 = Layer {
//             gates: vec![
//                 Gate {
//                     left: 0,
//                     right: 1,
//                     op: GateOp::Add,
//                     output: 0,
//                 },
//                 Gate {
//                     left: 2,
//                     right: 3,
//                     op: GateOp::Mul,
//                     output: 1,
//                 },
//             ],
//         };

//         let layer_3 = Layer {
//             gates: vec![Gate {
//                 left: 0,
//                 right: 1,
//                 op: GateOp::Add,
//                 output: 0,
//             }],
//         };

//         // let r_b = Fq::from(2);
//         // let r_c = Fq::from(3);

//         circuit.add_layer(layer_1);
//         circuit.add_layer(layer_2);
//         circuit.add_layer(layer_3);

//         let sub_claim = SubClaim {
//             challenges: vec![Fq::from(1), Fq::from(2)],
//             last_claimed_sum: Fq::from(0)
//         };

//         circuit.prove(sub_claim);
//     }
// }