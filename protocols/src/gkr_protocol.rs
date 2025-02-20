use std::vec;

use crate::{
    gkr_circuit::{Circuit, GateOp},
    multi_linear::MultiLinearPoly,
    partial_sum_check::{self, Proof},
    product_poly::ProductPoly,
    transcript::Transcript,
};
use ark_ff::PrimeField;

pub struct GKRProof<F: PrimeField> {
    pub output_layer: Vec<F>,    // an array of wᵢ
    pub w_i_evals: Vec<(F, F)>,  // array of wᵢ evaluated at r_b and r_c
    pub p_proofs: Vec<Proof<F>>, // array of sum-check proofs
}

impl<F: PrimeField> Circuit<F> {
    pub fn proof(&self) -> GKRProof<F> {
        let mut transcript = Transcript::new();
        let evaluated_circuit = self.evaluate();

        let mut sum_poly_array = Vec::new();
        let mut w_i_evals = Vec::new();
        let mut p_proofs = Vec::new();
        let mut r_a_challenges = Vec::new();

        let circuit_len = evaluated_circuit.len() - 1;

        // Get the output layer evaluations (W₀)
        let w_0 = evaluated_circuit[circuit_len].clone();

        // Pad W₀ to power of 2 if needed
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
        let output_layer = w_0_arr.clone();

        let w_0_len = w_0_arr.len().ilog2();

        // Get random point r₀
        transcript.absorb(&MultiLinearPoly::to_bytes(&w_0_arr));
        for _ in 0..w_0_len {
            let r_a = F::from_be_bytes_mod_order(&transcript.squeeze());
            r_a_challenges.push(r_a);
        }

        let w_0_eval = MultiLinearPoly::new(w_0_arr).evaluate(&r_a_challenges); // claimed sum = w_0(r)
        let init_claimed_sum = w_0_eval.computation[0];

        // f_ri_b_c = [add_i_ri_b_c * (w_i+1_b + w_i+1_c)] + [mul_i_ri_b_c * (w_i+1_b * w_i+1_c)]
        let next_layer_idx = circuit_len - 1;
        let (w_i_b_exploded, w_i_c_exploded) = self.explode_w_i(next_layer_idx);

        let sum_term = Circuit::<F>::element_wise_op(
            &w_i_b_exploded,
            &w_i_c_exploded,
            GateOp::Add,
        );
        let mul_term = Circuit::<F>::element_wise_op(&w_i_b_exploded, &w_i_c_exploded, GateOp::Mul);

        let (add_i, mul_i) = self.layer_i_add_mul(circuit_len);
        let mut add_i_mle = MultiLinearPoly::new(add_i);
        let mut mul_i_mle = MultiLinearPoly::new(mul_i);

        for r_a in r_a_challenges.iter() {
            add_i_mle = add_i_mle.partial_evaluate(*r_a, 0);
            mul_i_mle = mul_i_mle.partial_evaluate(*r_a, 0);
        }

        let p_poly_1 = ProductPoly::new(vec![
            add_i_mle,
            MultiLinearPoly {
                computation: sum_term,
            },
        ]);
        let p_poly_2 = ProductPoly::new(vec![
            mul_i_mle,
            MultiLinearPoly {
                computation: mul_term,
            },
        ]);

        let p_poly = vec![p_poly_1, p_poly_2];
        sum_poly_array.push(p_poly.clone());

        let p_proof = partial_sum_check::proof::<F>(p_poly, init_claimed_sum);
        p_proofs.push(p_proof.clone());
        let mut challenges = p_proof.challenges.clone();

        // For each layer i (going backwards from output to input)
        // since last layer has been done, we start with next layer
        // [0, 1, 2, 3] => would start at 2 and end at 1 as w will go down to 0
        for layer_idx in (1..circuit_len).rev() {
            let next_layer_idx = layer_idx - 1; // this is because w is 1 layer ahead
            let current_layer_w = evaluated_circuit[layer_idx].clone();

            // claimed_sum = (alpha * Wᵢ(*b)) + (beta * Wᵢ(*c))
            let claimed_sum = self.new_claimed_sum(current_layer_w, &challenges);

            // Get the add and mul vectors for current layer
            let (new_add, new_mul) = self.gkr_trick(&challenges, layer_idx);

            // Get the next layer evaluations (Wᵢ₊₁)
            let (w_i_b_exploded, w_i_c_exploded) = self.explode_w_i(next_layer_idx);

            // Compute f_rᵢ(b, c) = addᵢ(rᵢ,b,c)(Wᵢ₊₁(b) + Wᵢ₊₁(c)) + mulᵢ(rᵢ,b,c)(Wᵢ₊₁(b) * Wᵢ₊₁(c))
            let sum_term = Circuit::<F>::element_wise_op(
                &w_i_b_exploded,
                &w_i_c_exploded,
                GateOp::Add,
            );
            let mul_term =
                Circuit::<F>::element_wise_op(&w_i_b_exploded, &w_i_c_exploded, GateOp::Mul);

            // Create the polynomials for sum-check
            let p_poly_1 = ProductPoly::new(vec![
                new_add,
                MultiLinearPoly {
                    computation: sum_term,
                },
            ]);
            let p_poly_2 = ProductPoly::new(vec![
                new_mul,
                MultiLinearPoly {
                    computation: mul_term,
                },
            ]);

            let p_poly = vec![p_poly_1, p_poly_2];
            sum_poly_array.push(p_poly.clone());

            // Run sum-check protocol
            let p_proof = partial_sum_check::proof::<F>(p_poly, claimed_sum);
            p_proofs.push(p_proof.clone());

            challenges = p_proof.challenges.clone();
        }

        // this section is to get the evaluations of wᵢ at r_b and r_c to be used by the verifier
        for layer_idx in (0..circuit_len).rev() {
            let current_layer_w = evaluated_circuit[layer_idx].clone();
            let challenges = p_proofs[circuit_len - layer_idx - 1].challenges.clone();

            let mid = challenges.len() / 2;
            let (r_b_challenges, r_c_challenges) = challenges.split_at(mid);

            let w_i_b =
                MultiLinearPoly::new(current_layer_w.clone()).evaluate(&r_b_challenges);
            let w_i_c =
                MultiLinearPoly::new(current_layer_w.clone()).evaluate(&r_c_challenges);

            w_i_evals.push((w_i_b.computation[0], w_i_c.computation[0]));
        }

        GKRProof {
            output_layer,
            w_i_evals,
            p_proofs,
        }
    }

    pub fn verify(&self, proof: GKRProof<F>) -> bool {
        // recall that f(a, b, c) has already been evaluated by r_a to get f(b, c)
        // NOTE that the prover called evaluate meaning he has the Wᵢ values for every step while the verifier only has the input
        let mut transcript = Transcript::new();
        let mut last_challenges = Vec::new();
        let mut curr_challenges = Vec::new();
        let mut current_claimed_sum = F::zero();
        let circuit_len = self.layers.len(); // actual number of layers
        let mut last_idx = 0;

        let w_0_arr = proof.output_layer.clone();
        transcript.absorb(&MultiLinearPoly::to_bytes(&w_0_arr));
        let r_a = F::from_be_bytes_mod_order(&transcript.squeeze());

        let (add_i, mul_i) = self.layer_i_add_mul(circuit_len);
        let mut new_add = MultiLinearPoly::new(add_i).partial_evaluate(r_a, 0);
        let mut new_mul = MultiLinearPoly::new(mul_i).partial_evaluate(r_a, 0);

        for (i, p_proof) in proof.p_proofs.iter().enumerate() {
            let sub_claim = partial_sum_check::verify(p_proof.clone());
            let challenges = sub_claim.challenges.clone();

            curr_challenges = challenges.clone();
            last_idx = i;

            // For all but the last proof, check against w_i_evals
            if i < proof.p_proofs.len() - 1 {
                let new_add_eval = new_add.evaluate(&challenges);
                let new_mul_eval = new_mul.evaluate(&challenges);

                let (w_i_rb, w_i_rc) = proof.w_i_evals[i];
                let w_sum = w_i_rb + w_i_rc;
                let w_mul = w_i_rb * w_i_rc;

                let check =
                    (new_add_eval.computation[0] * w_sum) + (new_mul_eval.computation[0] * w_mul);

                if check != sub_claim.last_claimed_sum {
                    return false;
                }

                (new_add, new_mul) = self.gkr_trick(&challenges, circuit_len - i - 1);

                last_challenges = challenges.clone();
            }

            current_claimed_sum = sub_claim.last_claimed_sum;
        }

        // Finally, performs oracle check for each layer using the below
        // f(b, c) = [add_i(b, c) * (w_i+1(b) + w_i+1(c))] + [mul_i(b,c) * (w_i+1(b) * w_i+1(c))]
        let input_evaluations = self.inputs.clone();
        let mut input_poly = MultiLinearPoly::new(input_evaluations);

        let mid = curr_challenges.len() / 2;
        let (r_b_challenges, r_c_challenges) = curr_challenges.split_at(mid);

        let input_eval_b = input_poly.evaluate(&r_b_challenges).computation[0];
        let input_eval_c = input_poly.evaluate(&r_c_challenges).computation[0];
        let input_w_sum = input_eval_b + input_eval_c;
        let input_w_mul = input_eval_b * input_eval_c;

        (new_add, new_mul) = self.gkr_trick(&last_challenges, circuit_len - last_idx);
        let new_add_eval = new_add.evaluate(&curr_challenges).computation[0];
        let new_mul_eval = new_mul.evaluate(&curr_challenges).computation[0];

        let oracle_check = (new_add_eval * input_w_sum) + (new_mul_eval * input_w_mul);

        oracle_check == current_claimed_sum
    }
}

#[cfg(test)]
mod test {
    use crate::gkr_circuit::test::setup_test_circuit8;

    #[test]
    fn test_gkr_protocol_proof() {
        let circuit = setup_test_circuit8();

        circuit.proof();
    }

    #[test]
    fn test_gkr_protocol_verify() {
        let circuit = setup_test_circuit8();

        let proof = circuit.proof();
        let result = circuit.verify(proof);
        assert!(&result);
    }
}
