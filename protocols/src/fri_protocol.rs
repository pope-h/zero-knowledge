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
    pub values_at_index_proof: Vec<MerkleProof>,
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

    // -> Commitment<F>
    // This fn can be made to take in num_rounds in future impl
    pub fn generate_proof(&self) {
        let mut transcript = Transcript::new();
        let mut m_hashes = vec![];

        let mut f_poly = self.poly.clone();

        let padded_poly = self.pad_to_power_of_two();
        let fft = FastFourierTransform::new(padded_poly);
        let mut eval_poly = fft.evaluate().coefficients;

        let num_rounds = eval_poly.len().ilog2();

        for _i in 0..num_rounds {
            let poly_string: Vec<String> = eval_poly.iter().map(|d| d.to_string()).collect();
            let poly_bytes: Vec<&[u8]> = poly_string.iter().map(|s| s.as_bytes()).collect();
        
            let m_tree = MerkleTree::new(&poly_bytes);
            let m_root = m_tree.root().unwrap();
        
            transcript.absorb(&m_root);
            m_hashes.push(m_root);
        
            let r = F::from_be_bytes_mod_order(&transcript.squeeze());

            f_poly = fold_poly(&f_poly, r);
            let padded_poly = pad_poly_to_power_of_two(&f_poly);
        
            let fft = FastFourierTransform::new(padded_poly);
            eval_poly = fft.evaluate().coefficients;
        }

        //=========================================================================================
        // Sample a random index and get the evaluations at that index
        //=========================================================================================
        let verifier_field = F::from_be_bytes_mod_order(&transcript.squeeze());
        let field_integer_repr = verifier_field.into_bigint().as_ref()[0];
        let v_index = (field_integer_repr as usize) % self.poly.len();
        dbg!(&v_index);

        let values_at_index = self.poly[v_index];

        // Commitment { root_hashes, final_poly }
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