use ark_ff::{FftField, PrimeField};

use crate::{
    fri::fft::FastFourierTransform,
    fri::fri_helper_functions::{fold_poly, pad_poly_to_power_of_two},
    merkle_tree::{MerkleProof, MerkleTree},
    transcript::Transcript,
};

pub struct FRIProtocol<F: FftField> {
    pub poly: Vec<F>,
    pub blowup_factor: usize,
    pub max_degree: usize,
}

pub struct FRIProof<F: FftField> {
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
    pub fn new(poly: Vec<F>, blowup_factor: usize) -> Self {
        let max_degree = poly.len() - 1;
        FRIProtocol {
            poly,
            blowup_factor,
            max_degree,
        }
    }

    // This fn can be made to take in num_rounds in future impl
    pub fn generate_proof(&self) -> FRIProof<F> {
        let mut transcript = Transcript::new();
        let mut m_hashes = vec![];
        let mut m_trees = vec![];
        let mut c_sums = vec![];
        let mut v_at_index = vec![];
        let mut v_at_neg_index = vec![];
        let mut p_at_index = vec![];
        let mut p_at_neg_index = vec![];
        let mut all_evals = vec![];

        let mut f_poly = self.poly.clone();

        let padded_poly = self.pad_to_power_of_two();
        let domain_size = padded_poly.len();
        let fft = FastFourierTransform::new(padded_poly);
        let mut eval_poly = fft.evaluate().coefficients;

        all_evals.push(eval_poly.clone());

        let num_rounds = eval_poly.len().ilog2();

        for _i in 0..num_rounds {
            let poly_string: Vec<String> = eval_poly.iter().map(|d| d.to_string()).collect();
            let poly_bytes: Vec<&[u8]> = poly_string.iter().map(|s| s.as_bytes()).collect();

            let m_tree = MerkleTree::new(&poly_bytes);
            let m_root = m_tree.root().unwrap();

            transcript.absorb(&m_root);
            m_hashes.push(m_root);

            m_trees.push(m_tree);

            let r = F::from_be_bytes_mod_order(&transcript.squeeze());

            if f_poly.len() == 1 {
                eval_poly = fold_poly(&f_poly, r);

                all_evals.push(eval_poly.clone());

                break;
            } else {
                f_poly = fold_poly(&f_poly, r);
                let padded_poly = pad_poly_to_power_of_two(&f_poly);

                let fft = FastFourierTransform::new(padded_poly);
                eval_poly = fft.evaluate().coefficients;

                all_evals.push(eval_poly.clone());
            }
        }

        let final_poly = eval_poly;

        //=========================================================================================
        // Sample a random index and get the evaluations at that index
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

        FRIProof {
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

    pub fn verify(&self, proof: FRIProof<F>) -> bool {
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
        fri.generate_proof();
    }

    #[test]
    fn test_fri_protocol_verify() {
        let poly: Vec<ark_ff::Fp<ark_ff::MontBackend<ark_bn254::FrConfig, 4>, 4>> =
            vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let fri = FRIProtocol::new(poly, 2);
        let proof = fri.generate_proof();
        assert!(fri.verify(proof));
    }
}
