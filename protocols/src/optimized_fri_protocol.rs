use ark_ff::{FftField, PrimeField};

use crate::{fft::FastFourierTransform, fri_protocol::FRIProtocol, merkle_tree::MerkleTree, transcript::Transcript};

pub struct OptimizedFRIProof<F: FftField> {
    pub root_hashes: Vec<Vec<u8>>,
    pub final_poly: Vec<F>,
}

impl<F: FftField + PrimeField> FRIProtocol<F> {
    // -> Commitment<F>
    // This fn can be made to take in num_rounds in future impl
    pub fn optimized_proof(&self) {
        let mut transcript = Transcript::new();
        let mut m_hashes = vec![];

        let padded_poly = self.pad_to_power_of_two();
        let domain_size = padded_poly.len();

        //=========================================================================================
        // Get primitive root of unity for the domain
        //=========================================================================================
        let primitive_root = F::get_root_of_unity(domain_size as u64).unwrap();

        let fft = FastFourierTransform::new(padded_poly);
        let initial_evaluations = fft.evaluate().coefficients;

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

            transcript.absorb(&m_root);
            m_hashes.push(m_root);

            let r = F::from_be_bytes_mod_order(&transcript.squeeze());

            let next_domain_size = current_domain_size / 2;
            let mut next_evals = Vec::with_capacity(next_domain_size);

            for i in 0..next_domain_size {
                // Get the values at x and -x
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
            current_evals = next_evals;
            current_domain_size = next_domain_size;
            
            //=========================================================================================
            // Square the primitive root for the next domain
            // Since domain size is halved, the new primitive root is the square of the previous one
            //=========================================================================================
            current_primitive_root = current_primitive_root.square();
        }


        // // Return the final proof
        // FRIProof {
        //     root_hashes,
        //     final_poly: current_evals,
        // }
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
        fri.optimized_proof();
    }
}