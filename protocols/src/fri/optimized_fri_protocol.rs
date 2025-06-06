use std::vec;

use ark_ff::{FftField, PrimeField};

use crate::{
    fri::fft::FastFourierTransform,
    fri::fri_protocol::FRIProtocol,
    fri::merkle_tree::{MerkleProof, MerkleTree},
    transcript::Transcript,
};

pub struct OptimizedFRIProof<F: FftField> {
    pub root_hashes: Vec<Vec<u8>>,
    pub final_poly: Vec<F>,
    pub values_at_index: Vec<F>,
    pub values_at_neg_index: Vec<F>,
    pub merkle_trees: Vec<MerkleTree>,
    pub proofs_at_index: Vec<MerkleProof>,
    pub proofs_at_neg_index: Vec<MerkleProof>,
    pub claimed_sums: Vec<F>,
}

impl<F: FftField + PrimeField> FRIProtocol<F> {
    // This fn can be made to take in num_rounds in future impl
    pub fn optimized_proof(&self) -> OptimizedFRIProof<F> {
        let mut transcript = Transcript::new();
        let mut m_hashes = vec![];
        let mut m_trees = vec![];
        let mut c_sums = vec![];
        let mut v_at_index = vec![];
        let mut v_at_neg_index = vec![];
        let mut p_at_index = vec![];
        let mut p_at_neg_index = vec![];
        let mut all_evals = vec![];

        let padded_poly = self.pad_to_power_of_two();
        let domain_size = padded_poly.len();

        //=========================================================================================
        // Get primitive root of unity for the domain
        //=========================================================================================
        let primitive_root = F::get_root_of_unity(domain_size as u64).unwrap();

        let fft = FastFourierTransform::new(padded_poly);
        let initial_evaluations = fft.evaluate().coefficients;
        all_evals.push(initial_evaluations.clone());

        //=========================================================================================
        // Keep track of the current evaluations and domain parameters
        //=========================================================================================
        let mut current_evals = initial_evaluations;
        let mut current_domain_size = domain_size;
        let mut current_primitive_root = primitive_root;

        let num_rounds = domain_size.ilog2();

        for _round in 0..num_rounds {
            let poly_string: Vec<String> = current_evals.iter().map(|d| d.to_string()).collect();
            let poly_bytes: Vec<&[u8]> = poly_string.iter().map(|s| s.as_bytes()).collect();

            let m_tree = MerkleTree::new(&poly_bytes);
            let m_root = m_tree.root().unwrap();

            m_trees.push(m_tree);

            transcript.absorb(&m_root);
            m_hashes.push(m_root);

            let r = F::from_be_bytes_mod_order(&transcript.squeeze());

            let next_domain_size = current_domain_size / 2;
            let mut next_evals = Vec::with_capacity(next_domain_size);

            for i in 0..next_domain_size {
                //=========================================================================================
                // Get the values at x and -x
                //=========================================================================================
                let f_x = current_evals[i];
                let f_neg_x = current_evals[i + next_domain_size];

                //=========================================================================================
                // Get the actual domain element (ω^i)
                // i.e. root of unity raised to the power of i
                //=========================================================================================
                let omega_i = current_primitive_root.pow(&[i as u64]);

                //=========================================================================================
                // Calculate the next round value using the formula:
                // f₂(x²) = (f₁(x) + f₁(-x))/2 + r * ((f₁(x) - f₁(-x))/(2x))
                //=========================================================================================

                //=========================================================================================
                // First part: (f₁(x) + f₁(-x))/2
                //=========================================================================================
                let sum_term = (f_x + f_neg_x) * F::from(2).inverse().unwrap();

                //=========================================================================================
                // Second part: (f₁(x) - f₁(-x))/(2x)
                //=========================================================================================
                let diff = f_x - f_neg_x;
                let omega_i_doubled = omega_i.double();
                let diff_term = diff * omega_i_doubled.inverse().unwrap();

                //=========================================================================================
                // Final calculation
                //=========================================================================================
                let next_eval = sum_term + (r * diff_term);
                next_evals.push(next_eval);
            }

            //=========================================================================================
            // Update for next round
            //=========================================================================================
            all_evals.push(next_evals.clone());
            current_evals = next_evals;
            current_domain_size = next_domain_size;

            //=========================================================================================
            // Square the primitive root for the next domain
            // Since domain size is halved, the new primitive root is the square of the previous one
            //=========================================================================================
            current_primitive_root = current_primitive_root.square();
        }

        let final_poly = current_evals;

        //=========================================================================================
        // Sample a random index and get the evaluations at that index
        // This is the verifier's challenge
        //=========================================================================================
        let verifier_field = F::from_be_bytes_mod_order(&transcript.squeeze());
        let field_integer_repr = verifier_field.into_bigint().as_ref()[0];
        let mut v_index = (field_integer_repr as usize) % self.poly.len();
        let verifier_index = v_index;

        for round in 0..num_rounds {
            let round_domain_size = domain_size >> round;
            let half_domain_size = round_domain_size / 2;

            let value_at_index = all_evals[round as usize][v_index % round_domain_size];
            let proof_at_index =
                m_trees[round as usize].generate_proof(&value_at_index.to_string().as_bytes());
            v_at_index.push(value_at_index);
            p_at_index.push(proof_at_index.unwrap());

            let value_at_neg_index =
                all_evals[round as usize][(v_index + half_domain_size) % round_domain_size];
            let proof_at_neg_index =
                m_trees[round as usize].generate_proof(&value_at_neg_index.to_string().as_bytes());
            v_at_neg_index.push(value_at_neg_index);
            p_at_neg_index.push(proof_at_neg_index.unwrap());

            //=========================================================================================
            // We skip the first round since there is no claimed sum in it
            // The claimed_sum is computed so the verifier can have a direct comparison
            //=========================================================================================
            if round != 0 {
                let claimed_sum = all_evals[round as usize][verifier_index % round_domain_size];
                c_sums.push(claimed_sum);
            }

            v_index /= 2;
        }

        OptimizedFRIProof {
            root_hashes: m_hashes,
            final_poly,
            values_at_index: v_at_index,
            values_at_neg_index: v_at_neg_index,
            merkle_trees: m_trees,
            proofs_at_index: p_at_index,
            proofs_at_neg_index: p_at_neg_index,
            claimed_sums: c_sums,
        }
    }

    pub fn optimized_verify(&self, proof: OptimizedFRIProof<F>) -> bool {
        let mut transcript = Transcript::new();

        let root_hashes = proof.root_hashes;
        let values_at_index = proof.values_at_index;
        let values_at_neg_index = proof.values_at_neg_index;
        let merkle_trees = proof.merkle_trees;
        let proofs_at_index = proof.proofs_at_index;
        let proofs_at_neg_index = proof.proofs_at_neg_index;
        let claimed_sums = proof.claimed_sums;

        let domain_size = 2u64.pow(root_hashes.len() as u32);

        //=========================================================================================
        // Get primitive root of unity for the domain
        //=========================================================================================
        let mut primitive_root = F::get_root_of_unity(domain_size as u64).unwrap();

        let num_rounds = root_hashes.len();

        for index in 0..(num_rounds - 1) {
            let check_proof_i = merkle_trees[index].verify_proof(
                &values_at_index[index].to_string().as_bytes(),
                &proofs_at_index[index],
                &root_hashes[index],
            );
            let check_proof_neg_i = merkle_trees[index].verify_proof(
                &values_at_neg_index[index].to_string().as_bytes(),
                &proofs_at_neg_index[index],
                &root_hashes[index],
            );

            if !check_proof_i && !check_proof_neg_i {
                return false;
            }

            transcript.absorb(&root_hashes[index]);
            let r = F::from_be_bytes_mod_order(&transcript.squeeze());

            //=========================================================================================
            // Get the values at x and -x
            //=========================================================================================
            let f_x = values_at_index[index];
            let f_neg_x = values_at_neg_index[index];

            //=========================================================================================
            // Get the actual domain element (ω^i)
            // i.e. root of unity raised to the power of i
            //=========================================================================================
            let omega_i = primitive_root.pow(&[proofs_at_index[index].leaf_index as u64]);

            //=========================================================================================
            // Calculate the next round value using the formula:
            // f₂(x²) = (f₁(x) + f₁(-x))/2 + r * ((f₁(x) - f₁(-x))/(2x))
            //=========================================================================================

            //=========================================================================================
            // First part: (f₁(x) + f₁(-x))/2
            //=========================================================================================
            let sum_term = (f_x + f_neg_x) * F::from(2).inverse().unwrap();

            //=========================================================================================
            // Second part: (f₁(x) - f₁(-x))/(2x)
            //=========================================================================================
            let diff = f_x - f_neg_x;
            let omega_i_doubled = omega_i.double();
            let diff_term = diff * omega_i_doubled.inverse().unwrap();

            //=========================================================================================
            // Final calculation
            //=========================================================================================
            let expected_next_eval = sum_term + (r * diff_term);

            if claimed_sums[index] != expected_next_eval {
                return false;
            }

            primitive_root = primitive_root.square();
        }

        //=========================================================================================
        // Oracle check for the last round
        //=========================================================================================
        let check_proof_last_i = merkle_trees[num_rounds - 1].verify_proof(
            &values_at_index[num_rounds - 1].to_string().as_bytes(),
            &proofs_at_index[num_rounds - 1],
            &root_hashes[num_rounds - 1],
        );
        let check_proof_neg_last_i = merkle_trees[num_rounds - 1].verify_proof(
            &values_at_neg_index[num_rounds - 1].to_string().as_bytes(),
            &proofs_at_neg_index[num_rounds - 1],
            &root_hashes[num_rounds - 1],
        );

        if !check_proof_last_i && !check_proof_neg_last_i {
            return false;
        }

        transcript.absorb(&root_hashes[num_rounds - 1]);
        let r = F::from_be_bytes_mod_order(&transcript.squeeze());

        let f_x = values_at_index[num_rounds - 1];
        let f_neg_x = values_at_neg_index[num_rounds - 1];

        let omega_i = primitive_root.pow(&[proofs_at_index[num_rounds - 1].leaf_index as u64]);

        let sum_term = (f_x + f_neg_x) * F::from(2).inverse().unwrap();

        let diff = f_x - f_neg_x;
        let omega_i_doubled = omega_i.double();
        let diff_term = diff * omega_i_doubled.inverse().unwrap();

        let expected_last_eval = sum_term + (r * diff_term);

        proof.final_poly[0] == expected_last_eval
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr;

    #[test]
    fn test_fri_protocol() {
        let poly: Vec<ark_ff::Fp<ark_ff::MontBackend<ark_bn254::FrConfig, 4>, 4>> =
            vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let fri = FRIProtocol::new(poly, 2);
        fri.optimized_proof();
    }

    #[test]
    fn test_fri_protocol_verify() {
        let poly: Vec<ark_ff::Fp<ark_ff::MontBackend<ark_bn254::FrConfig, 4>, 4>> =
            vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let fri = FRIProtocol::new(poly, 2);
        let proof = fri.optimized_proof();

        assert!(fri.optimized_verify(proof));
    }
}
