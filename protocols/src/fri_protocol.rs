use ark_ff::{FftField, PrimeField};

use crate::{fft::FastFourierTransform, fri_helper_functions::fold_poly, merkle_tree::MerkleTree, transcript::Transcript};

pub struct FRIProtocol<F: FftField> {
    pub poly: Vec<F>,
    pub blowup_factor: usize,
    pub max_degree: usize,
}

pub struct FRIProof<F: FftField> {
    pub root_hashes: Vec<Vec<u8>>,
    pub final_poly: Vec<F>,
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
        // let this_poly = FRIProtocol::new(poly, blowup_factor);
        let mut transcript = Transcript::new();
        let mut m_hashes = vec![];

        let padded_poly = self.pad_to_power_of_two();
        dbg!(padded_poly.len());
        let fft = FastFourierTransform::new(padded_poly);
        let mut eval_poly = fft.evaluate().coefficients;
        dbg!(&eval_poly.len());

        let num_rounds = eval_poly.len().ilog2();
        dbg!(num_rounds);

        for _i in 0..num_rounds {
            dbg!(&eval_poly);
            let poly_string: Vec<String> = eval_poly.iter().map(|d| d.to_string()).collect();
            let poly_bytes: Vec<&[u8]> = poly_string.iter().map(|s| s.as_bytes()).collect();

            let m_tree = MerkleTree::new(&poly_bytes);
            let m_root = m_tree.root().unwrap();

            transcript.absorb(&m_root);
            m_hashes.push(m_root);

            dbg!("Got here");
            let r = F::from_be_bytes_mod_order(&transcript.squeeze());
            dbg!(&r);
            let current_poly = fold_poly(&eval_poly, r);

            let fft = FastFourierTransform::new(current_poly);
            eval_poly = fft.evaluate().coefficients;
            dbg!(&eval_poly.len());
        }

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