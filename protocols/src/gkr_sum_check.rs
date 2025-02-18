use ark_ff::PrimeField;
use crate::{gkr_protocol::{Circuit, GateOp}, multi_linear::MultiLinearPoly, partial_sum_check::{self}, product_poly::ProductPoly, transcript::Transcript};

pub struct GKRProof<F: PrimeField> {
    pub sum_poly_array: Vec<Vec<ProductPoly<F>>>,
    pub claimed_sums: Vec<F>,
}

impl<F: PrimeField> Circuit<F> {
    pub fn proof(&self) -> GKRProof<F> {
        let mut transcript = Transcript::new();
        let evaluated_circuit = self.evaluate();

        let mut sum_poly_array = Vec::new();
        let mut claimed_sums = Vec::new();
        
        // Get the output layer evaluations (W₀)
        let w_0 = evaluated_circuit[evaluated_circuit.len() - 1].clone();
        
        // Pad w_0 to power of 2 if needed
        let w_0_arr = if w_0.len() == 1 {
            vec![w_0[0], F::zero()]
        } else if w_0.len().is_power_of_two() {
            w_0
        } else {
            let target_length = w_0.len().next_power_of_two();
            let mut padded = w_0.clone();
            padded.resize(target_length, F::zero());
            padded
        };

        // Get random point r₀
        transcript.absorb(&MultiLinearPoly::to_bytes(w_0_arr.clone()));
        let r_a = F::from_be_bytes_mod_order(&transcript.squeeze());

        let w_0_eval = MultiLinearPoly::new(w_0_arr).partial_evaluate(r_a, 0); // claimed sum = w_0(r)
        claimed_sums.push(w_0_eval.computation[0]);

        // f_ri_b_c = [add_i_ri_b_c * (w_i+1_b + w_i+1_c)] + [mul_i_ri_b_c * (w_i+1_b * w_i+1_c)]
        let next_layer_w = evaluated_circuit.len() - 2;
        let (w_i_b_exploded, w_i_c_exploded) = self.eval_layer_i(next_layer_w);

        let sum_term = Circuit::<F>::element_wise_op(w_i_b_exploded.clone(), w_i_c_exploded.clone(), GateOp::Add);
        let mul_term = Circuit::<F>::element_wise_op(w_i_b_exploded, w_i_c_exploded, GateOp::Mul);

        let (add_i, mul_i) = self.layer_i_add_mul(0);
        let add_i_ri = MultiLinearPoly::new(add_i).partial_evaluate(r_a, 0);
        let mul_i_ri = MultiLinearPoly::new(mul_i).partial_evaluate(r_a, 0);

        let p_poly_1 = ProductPoly::new(vec![add_i_ri, MultiLinearPoly { computation: sum_term }]);
        let p_poly_2 = ProductPoly::new(vec![mul_i_ri, MultiLinearPoly { computation: mul_term }]);

        let p_poly = vec![p_poly_1, p_poly_2];
        sum_poly_array.push(p_poly.clone());

        let p_proof = partial_sum_check::proof::<F>(p_poly, w_0_eval.computation[0]);
        let challenges = p_proof.challenges;

        // For each layer i (going backwards from output to input)
        // since layer 0 has been done above, we start with layer 1
        // i am thinking 1..(self.layers.len() - 1) better so for e.g. 
        // [0, 1, 2, 3] => would start at 2 and end at 1 as w will go down to 0
        for layer_idx in (1..(self.layers.len() - 1)).rev() {
            let next_layer_w = layer_idx + 1;   // this is because w is 1 layer ahead
            let current_layer_w = evaluated_circuit[layer_idx].clone();

            // claimed_sum = (alpha * w_i(*b)) + (beta * w_i(*c))
            let claimed_sum = self.new_claimed_sum(current_layer_w, challenges.clone());
            claimed_sums.push(claimed_sum);

            // Get the add and mul vectors for current layer
            let (new_add, new_mul) = self.gkr_trick(challenges.clone(), layer_idx);

            // Evaluate Wi+1 at points b* and c*
            let (w_i_b_exploded, w_i_c_exploded) = self.eval_layer_i(next_layer_w);

            // Compute f_ri(b, c) = add_i(ri,b,c)(Wi+1(b) + Wi+1(c)) + mul_i(ri,b,c)(Wi+1(b) * Wi+1(c))
            let sum_term = Circuit::<F>::element_wise_op(
                w_i_b_exploded.clone(), 
                w_i_c_exploded.clone(), 
                GateOp::Add
            );
            let mul_term = Circuit::<F>::element_wise_op(
                w_i_b_exploded, 
                w_i_c_exploded, 
                GateOp::Mul
            );

            // Create the polynomials for sum-check
            let p_poly_1 = ProductPoly::new(vec![new_add, MultiLinearPoly { computation: sum_term }]);
            let p_poly_2 = ProductPoly::new(vec![new_mul, MultiLinearPoly { computation: mul_term }]);

            let p_poly = vec![p_poly_1, p_poly_2];
            sum_poly_array.push(p_poly.clone());

            // Run sum-check protocol
            let p_proof = partial_sum_check::proof::<F>(p_poly, claimed_sum);
            dbg!(&p_proof);

            // Verify sum-check protocol. Would be removed, jus for sanity check
            let p_verify = partial_sum_check::verify(p_proof);
            dbg!(&p_verify);
        }

        GKRProof {
            sum_poly_array,
            claimed_sums,
        }
    }
}

#[cfg(test)]
mod test {
    use ark_bn254::Fq;
    use crate::gkr_protocol::{Circuit, Gate, GateOp, Layer};

    #[test]
    fn setup_test_circuit() {
        let inputs = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
        ];
        let mut circuit = Circuit::new(inputs);

        let layer_1 = Layer {
            gates: vec![
                Gate {
                    left: 0,
                    right: 1,
                    op: GateOp::Add,
                    output: 0,
                },
                Gate {
                    left: 2,
                    right: 3,
                    op: GateOp::Mul,
                    output: 1,
                },
                Gate {
                    left: 4,
                    right: 5,
                    op: GateOp::Mul,
                    output: 2,
                },
                Gate {
                    left: 6,
                    right: 7,
                    op: GateOp::Mul,
                    output: 3,
                },
            ],
        };

        let layer_2 = Layer {
            gates: vec![
                Gate {
                    left: 0,
                    right: 1,
                    op: GateOp::Add,
                    output: 0,
                },
                Gate {
                    left: 2,
                    right: 3,
                    op: GateOp::Mul,
                    output: 1,
                },
            ],
        };

        let layer_3 = Layer {
            gates: vec![Gate {
                left: 0,
                right: 1,
                op: GateOp::Add,
                output: 0,
            }],
        };

        circuit.add_layer(layer_1);
        circuit.add_layer(layer_2);
        circuit.add_layer(layer_3);

        circuit.proof();
    }
}