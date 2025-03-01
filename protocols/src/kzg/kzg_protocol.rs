use ark_ec::{pairing::Pairing, PrimeGroup};
use ark_ff::{PrimeField, Zero};

use crate::{kzg::kzg_helper_functions::{blow_up, compute_commitment, compute_poly_minus_v, compute_quotient, compute_remainder}, multi_linear::MultiLinearPoly};

pub struct KZGProof<P: Pairing> {
    pub quotient_evals: Vec<P::G1>,
}

pub fn proof<F: PrimeField, P: Pairing>(poly: MultiLinearPoly<F>, encrypted_basis: &[P::G1], vars_to_open: &[F]) -> KZGProof<P> {
    // Since we can't evaluate Q(x) at tau as we don't know tau, 
    // we need to first blow up Q_a(b, c) or Q_b(c) back to Q(a, b, c)
    // Then we mul and add of Q(a, b, c) with the lagrange_basis

    let mut quotient_evals = Vec::new();
    // let v = poly.evaluate(vars_to_open).computation[0];

    let commitment = compute_commitment::<F, P>(&poly, encrypted_basis);
    dbg!(&commitment);

    let mut poly_minus_v = compute_poly_minus_v(poly, vars_to_open);

    for i in 0..(vars_to_open.len()) {
        let quotient = compute_quotient(&poly_minus_v);
        let blown_quotient = blow_up(quotient, i + 1);

        let mut quotient_eval = P::G1::zero();
        for (j, e_basis) in encrypted_basis.iter().enumerate() {
            quotient_eval += e_basis.mul_bigint(blown_quotient.computation[j].into_bigint());
        }
        quotient_evals.push(quotient_eval);

        let remainder = compute_remainder(poly_minus_v, vars_to_open[i]);
        poly_minus_v = remainder;
    }
    assert_eq!(poly_minus_v.computation[0], F::zero());

    KZGProof { quotient_evals }
}

#[cfg(test)]
mod test {
    use crate::{kzg::kzg_helper_functions::test::poly_1, kzg::trusted_setup::tests::setup};

    use super::*;
    use ark_bls12_381::{Bls12_381, Fr as BlsFr};

    #[test]
    fn test_proof() {
        let setup = setup();
        let poly = poly_1();
        let vars_to_open = vec![BlsFr::from(6), BlsFr::from(4), BlsFr::from(0)];

        let proof = proof::<BlsFr, Bls12_381>(poly, &setup.g1_arr, &vars_to_open);
        assert_eq!(proof.quotient_evals.len(), 3);
    }
}