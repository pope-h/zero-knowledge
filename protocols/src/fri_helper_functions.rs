use ark_ff::FftField;

use crate::fri_protocol::FRIProtocol;

impl<F: FftField> FRIProtocol<F> {
    pub fn pad_to_power_of_two(&self) -> Vec<F> {
        let mut padded_poly = self.poly.clone();
        let domain_size = self.domain_size();
        padded_poly.resize(domain_size, F::zero());

        padded_poly
    }

    //=========================================================================================
    // This functions ensure that the domain size is a power of 2
    // i.e. if num is 3 with blow_up_factor 2, the domain size will be 8 and not 6
    //=========================================================================================
    pub fn domain_size(&self) -> usize {
        let min_size = (self.max_degree + 1) * self.blowup_factor;
        let mut size = 1;

        while size < min_size {
            size *= 2;
        }

        size
    }

    pub fn fold_poly(&self, r_challenge: F) -> Vec<F> {
        let (even, odd) = self.split_poly();

        even.iter()
            .zip(odd.iter())
            .map(|(e, o)| *e + (r_challenge * *o))
            .collect()
    }

    fn split_poly(&self) -> (Vec<F>, Vec<F>) {
        let mut even = Vec::new();
        let mut odd = Vec::new();

        for (i, coeff) in self.poly.iter().enumerate() {
            if i % 2 == 0 {
                even.push(*coeff);
            } else {
                odd.push(*coeff);
            }
        }

        (even, odd)
    }

    pub fn pad_vector_to_power_of_two(&self, vec: Vec<F>) -> Vec<F> {
        let mut padded_vec = vec;
        let domain_size = self.domain_size();
        padded_vec.resize(domain_size, F::zero());

        padded_vec
    }
}

#[cfg(test)]
mod test {
    use crate::fri_protocol::FRIProtocol;
    use ark_bn254::Fq;

    fn poly_1() -> Vec<Fq> {
        vec![Fq::from(1), Fq::from(2), Fq::from(3)]
    }

    #[test]
    fn test_pad_to_power_of_two() {
        let poly = poly_1();
        let fri_protocol = FRIProtocol::new(poly, 2);

        let padded_poly = fri_protocol.pad_to_power_of_two();
        assert_eq!(padded_poly.len(), 8);
    }

    #[test]
    fn test_domain_size() {
        let poly = poly_1();
        let fri_protocol = FRIProtocol::new(poly, 2);

        let domain_size = fri_protocol.domain_size();
        assert_eq!(domain_size, 8);
    }

    #[test]
    fn test_fold_poly() {
        let poly = poly_1();
        let fri_protocol = FRIProtocol::new(poly, 2);

        let r_challenge = Fq::from(2);
        let folded_poly = fri_protocol.fold_poly(r_challenge);

        // (0 + 2*2) + (1 + 2*0) => 5
        let expected_folded_poly = vec![Fq::from(5)];

        assert_eq!(folded_poly, expected_folded_poly);
    }

    #[test]
    fn test_pad_vector_to_power_of_two() {
        let poly = poly_1();
        let fri_protocol = FRIProtocol::new(poly, 2);

        let vec = vec![Fq::from(5), Fq::from(6), Fq::from(7)];
        let padded_vec = fri_protocol.pad_vector_to_power_of_two(vec);

        let expected_padded_vec = vec![
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
        ];

        assert_eq!(padded_vec, expected_padded_vec);
    }
}