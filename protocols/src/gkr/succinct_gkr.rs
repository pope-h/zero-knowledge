use crate::{
    gkr::{gkr_circuit::Circuit, partial_sum_check::Proof},
    kzg::kzg_helper_functions::{
        blow_up, compute_commitment, compute_poly_minus_v, compute_quotient, compute_remainder,
    },
    multi_linear::MultiLinearPoly,
    transcript::Transcript,
};
use ark_ec::{
    pairing::{Pairing, PairingOutput},
    PrimeGroup,
};
use ark_ff::{AdditiveGroup, PrimeField, Zero};

use super::{gkr_circuit::GateOp, partial_sum_check, product_poly::ProductPoly};

#[derive(Debug)]
pub struct SuccinctGKRProof<F: PrimeField, P: Pairing> {
    pub output_layer: Vec<F>,    // an array of wᵢ
    pub w_i_evals: Vec<(F, F)>,  // array of wᵢ evaluated at r_b and r_c
    pub p_proofs: Vec<Proof<F>>, // array of sum-check proofs
    pub commitment: P::G1,
    pub quotient_evals_rb: Vec<P::G1>,
    pub quotient_evals_rc: Vec<P::G1>,
}

impl<F: PrimeField> Circuit<F> {
    pub fn succinct_proof<P: Pairing>(&self, encrypted_basis: &[P::G1]) -> SuccinctGKRProof<F, P> {
        let mut transcript = Transcript::new();
        let evaluated_circuit = self.evaluate();

        let mut sum_poly_array = Vec::new();
        let mut w_i_evals = Vec::new();
        let mut p_proofs = Vec::new();
        let mut r_a_challenges = Vec::new();

        //=========================================================================================
        // GKR Proving Process
        //=========================================================================================
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

        let w_0_eval = MultiLinearPoly::new(&w_0_arr).evaluate(&r_a_challenges); // claimed sum = w_0(r)
        let init_claimed_sum = w_0_eval.computation[0];

        //=========================================================================================
        // f_ri_b_c = [add_i_ri_b_c * (w_i+1_b + w_i+1_c)] + [mul_i_ri_b_c * (w_i+1_b * w_i+1_c)]
        //=========================================================================================
        let next_layer_idx = circuit_len - 1;
        let (w_i_b_exploded, w_i_c_exploded) = self.explode_w_i(next_layer_idx);

        let sum_term = Circuit::<F>::element_wise_op(&w_i_b_exploded, &w_i_c_exploded, GateOp::Add);
        let mul_term = Circuit::<F>::element_wise_op(&w_i_b_exploded, &w_i_c_exploded, GateOp::Mul);

        let (add_i, mul_i) = self.layer_i_add_mul(circuit_len);
        let mut add_i_mle = MultiLinearPoly::new(&add_i);
        let mut mul_i_mle = MultiLinearPoly::new(&mul_i);

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

        //=========================================================================================
        // For each layer i (going backwards from output to input)
        // since last layer has been done, we start with next layer
        // [0, 1, 2, 3] => would start at 2 and end at 1 as w will go down to 0
        //=========================================================================================
        for layer_idx in (1..circuit_len).rev() {
            let next_layer_idx = layer_idx - 1; // this is because w is 1 layer ahead
            let current_layer_w = evaluated_circuit[layer_idx].clone();

            // claimed_sum = (alpha * Wᵢ(*b)) + (beta * Wᵢ(*c))
            let claimed_sum = self.new_claimed_sum(current_layer_w, &challenges);

            // Get the add and mul vectors for current layer
            let (new_add, new_mul) = self.gkr_trick(&challenges, layer_idx);

            // Get the next layer evaluations (Wᵢ₊₁)
            let (w_i_b_exploded, w_i_c_exploded) = self.explode_w_i(next_layer_idx);

            //=========================================================================================
            // Compute f_rᵢ(b, c) = addᵢ(rᵢ,b,c)(Wᵢ₊₁(b) + Wᵢ₊₁(c)) + mulᵢ(rᵢ,b,c)(Wᵢ₊₁(b) * Wᵢ₊₁(c))
            //=========================================================================================
            let sum_term =
                Circuit::<F>::element_wise_op(&w_i_b_exploded, &w_i_c_exploded, GateOp::Add);
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

        //=========================================================================================
        // GKR evaluations of wᵢ at r_b and r_c to be used by the verifier
        //=========================================================================================
        for layer_idx in (0..circuit_len).rev() {
            let current_layer_w = evaluated_circuit[layer_idx].clone();
            let challenges = p_proofs[circuit_len - layer_idx - 1].challenges.clone();

            let mid = challenges.len() / 2;
            let (r_b_challenges, r_c_challenges) = challenges.split_at(mid);

            let w_i_b = MultiLinearPoly::new(&current_layer_w).evaluate(&r_b_challenges);
            let w_i_c = MultiLinearPoly::new(&current_layer_w).evaluate(&r_c_challenges);

            // transcript.absorb(&MultiLinearPoly::to_bytes(&[w_i_b, w_i_c]));

            w_i_evals.push((w_i_b.computation[0], w_i_c.computation[0]));
        }

        //=========================================================================================
        // KZG Proof
        //=========================================================================================
        let mut quotient_evals_rb = Vec::new();
        let mut quotient_evals_rc = Vec::new();
        let input_poly = MultiLinearPoly::new(&self.inputs);

        let commitment = compute_commitment::<F, P>(&input_poly, encrypted_basis);

        // Get last challenges for the quotient evaluations
        let final_challenges = p_proofs.last().unwrap().challenges.clone();
        let mid = final_challenges.len() / 2;
        let (r_b_challenges, r_c_challenges) = final_challenges.split_at(mid);

        //=========================================================================================
        // Generating Quotients Q(τ) for r_b
        //=========================================================================================
        let mut poly_minus_v_b = compute_poly_minus_v(input_poly.clone(), &r_b_challenges);
        for i in 0..(r_b_challenges.len()) {
            let quotient = compute_quotient(&poly_minus_v_b);
            let blown_quotient = blow_up(quotient, i + 1);

            let mut quotient_eval = P::G1::zero();
            for (j, e_basis) in encrypted_basis.iter().enumerate() {
                quotient_eval += e_basis.mul_bigint(blown_quotient.computation[j].into_bigint());
            }
            quotient_evals_rb.push(quotient_eval);

            let remainder = compute_remainder(poly_minus_v_b, r_b_challenges[i]);
            poly_minus_v_b = remainder;
        }

        //=========================================================================================
        // Generating Quotients Q(τ) for r_c
        //=========================================================================================
        let mut poly_minus_v_c = compute_poly_minus_v(input_poly, &r_c_challenges);
        for i in 0..(r_c_challenges.len()) {
            let quotient = compute_quotient(&poly_minus_v_c);
            let blown_quotient = blow_up(quotient, i + 1);

            let mut quotient_eval = P::G1::zero();
            for (j, e_basis) in encrypted_basis.iter().enumerate() {
                quotient_eval += e_basis.mul_bigint(blown_quotient.computation[j].into_bigint());
            }
            quotient_evals_rc.push(quotient_eval);

            let remainder = compute_remainder(poly_minus_v_c, r_c_challenges[i]);
            poly_minus_v_c = remainder;
        }

        SuccinctGKRProof {
            output_layer,
            w_i_evals,
            p_proofs,
            commitment,
            quotient_evals_rb,
            quotient_evals_rc,
        }
    }

    pub fn succinct_verify<P: Pairing>(
        &self,
        proof: &SuccinctGKRProof<F, P>,
        encrypted_basis_g2: &[P::G2],
    ) -> bool {
        let g1_generator = P::G1::generator();
        let g2_generator = P::G2::generator();

        let mut transcript = Transcript::new();
        let mut last_challenges = Vec::new();
        let mut curr_challenges = Vec::new();
        let mut current_claimed_sum = F::zero();
        let circuit_len = self.layers.len(); // actual number of layers
        let mut last_idx = 0;

        //=========================================================================================
        // GKR Verification Process
        //=========================================================================================
        let w_0_arr = proof.output_layer.clone();
        transcript.absorb(&MultiLinearPoly::to_bytes(&w_0_arr));
        let r_a = F::from_be_bytes_mod_order(&transcript.squeeze());

        let (add_i, mul_i) = self.layer_i_add_mul(circuit_len);
        let mut new_add = MultiLinearPoly::new(&add_i).partial_evaluate(r_a, 0);
        let mut new_mul = MultiLinearPoly::new(&mul_i).partial_evaluate(r_a, 0);

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

        let mid = curr_challenges.len() / 2;
        let (r_b_challenges, r_c_challenges) = curr_challenges.split_at(mid);

        let input_eval_b = proof.w_i_evals.last().unwrap().0;
        let input_eval_c = proof.w_i_evals.last().unwrap().1;

        //=========================================================================================
        // KZG Verification Process
        // Verify opening at r_b
        // pairing(g1_(f(τ) - v), g2_1) == pairing(Σ(g1_Q(τ), g2_(τ - a)))
        //=========================================================================================
        let b_lhs = P::pairing(
            proof.commitment - g1_generator.mul_bigint(input_eval_b.into_bigint()),
            g2_generator.mul_bigint(F::one().into_bigint()),
        );

        let mut b_rhs = PairingOutput::ZERO;
        for (i, tau) in encrypted_basis_g2.iter().enumerate() {
            b_rhs += P::pairing(
                proof.quotient_evals_rb[i],
                *tau - g2_generator.mul_bigint(r_b_challenges[i].into_bigint()),
            );
        }

        // r_b check
        if b_lhs != b_rhs {
            return false;
        }

        //=========================================================================================
        // KZG Verification Process
        // Verify opening at r_c
        // pairing(g1_(f(τ) - v), g2_1) == pairing(Σ(g1_Q(τ), g2_(τ - a)))
        //=========================================================================================
        let c_lhs = P::pairing(
            proof.commitment - g1_generator.mul_bigint(input_eval_c.into_bigint()),
            g2_generator.mul_bigint(F::one().into_bigint()),
        );

        let mut c_rhs = PairingOutput::ZERO;
        for (i, tau) in encrypted_basis_g2.iter().enumerate() {
            c_rhs += P::pairing(
                proof.quotient_evals_rc[i],
                *tau - g2_generator.mul_bigint(r_c_challenges[i].into_bigint()),
            );
        }

        // r_c check
        if c_lhs != c_rhs {
            return false;
        }

        //=========================================================================================
        // Input layer is verified, now perform the GKR oracle check
        // f(b, c) = [add_i(b, c) * (w_i+1(b) + w_i+1(c))] + [mul_i(b,c) * (w_i+1(b) * w_i+1(c))]
        //=========================================================================================
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
    use crate::{
        gkr::gkr_circuit::{Circuit, Gate, GateOp, Layer},
        kzg::trusted_setup::tests::setup,
    };

    use ark_bls12_381::{Bls12_381, Fr as BlsFr};

    pub fn setup_test_circuit_s() -> Circuit<BlsFr> {
        let inputs = vec![
            BlsFr::from(1),
            BlsFr::from(2),
            BlsFr::from(3),
            BlsFr::from(4),
            BlsFr::from(5),
            BlsFr::from(6),
            BlsFr::from(7),
            BlsFr::from(8),
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

        circuit
    }

    #[test]
    fn test_gkr_protocol_proof() {
        let circuit = setup_test_circuit_s();
        let setup = setup();

        // let result = circuit.succinct_proof::<Bls12_381>(&setup.g1_arr);
        // dbg!(&result);

        circuit.succinct_proof::<Bls12_381>(&setup.g1_arr);
    }

    #[test]
    fn test_gkr_protocol_verify() {
        let circuit = setup_test_circuit_s();
        let setup = setup();

        let proof = circuit.succinct_proof::<Bls12_381>(&setup.g1_arr);
        let result = circuit.succinct_verify::<Bls12_381>(&proof, &setup.g2_arr);

        assert!(&result);
    }
}
