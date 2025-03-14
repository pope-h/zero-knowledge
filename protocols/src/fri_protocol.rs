use ark_ff::{FftField, PrimeField};

use crate::{fft::FastFourierTransform, fri_helper_functions::{fold_poly, pad_poly_to_power_of_two}, merkle_tree::{MerkleProof, MerkleTree}, transcript::Transcript};

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
    pub proofs_at_index: Vec<MerkleProof>,
    pub proofs_at_neg_index: Vec<MerkleProof>,
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

        for round in 0..num_rounds{
            let round_domain_size = domain_size >> round;
            let half_domain_size = round_domain_size / 2;

            let value_at_index = all_evals[round as usize][v_index % round_domain_size];
            let proof_at_index = m_trees[round as usize].generate_proof(&value_at_index.to_string().as_bytes());
            v_at_index.push(value_at_index);
            p_at_index.push(proof_at_index.unwrap());

            let value_at_neg_index = all_evals[round as usize][(v_index + half_domain_size) % round_domain_size];
            let proof_at_neg_index = m_trees[round as usize].generate_proof(&value_at_neg_index.to_string().as_bytes());
            v_at_neg_index.push(value_at_neg_index);
            p_at_neg_index.push(proof_at_neg_index.unwrap());

            v_index /= 2;
        }

        FRIProof {
            root_hashes: m_hashes,
            final_poly,
            values_at_index: v_at_index,
            values_at_neg_index: v_at_neg_index,
            proofs_at_index: p_at_index,
            proofs_at_neg_index: p_at_neg_index,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr;

    #[test]
    fn test_fri_protocol() {
        let poly: Vec<ark_ff::Fp<ark_ff::MontBackend<ark_bn254::FrConfig, 4>, 4>> = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let fri = FRIProtocol::new(poly, 2);
        fri.generate_proof();
    }
}