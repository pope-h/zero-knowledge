use ark_ff::FftField;

// use crate::transcript::Transcript;

pub struct FRIProtocol<F: FftField> {
    pub poly: Vec<F>,
    pub blowup_factor: usize,
    pub max_degree: usize,
}

pub struct FRIProof<F: FftField> {
    pub root_hashes: Vec<Vec<u8>>,
    pub final_poly: Vec<F>,
}

impl<F: FftField> FRIProtocol<F> {
    pub fn new(poly: Vec<F>, blowup_factor: usize) -> Self {
        let max_degree = poly.len() - 1;
        FRIProtocol {
            poly,
            blowup_factor,
            max_degree,
        }
    }

    // pub fn generate_proof(&self) -> Commitment<F> {
    //     let transcript = Transcript::new();

    //     Commitment { root_hashes, final_poly }
    // }
}