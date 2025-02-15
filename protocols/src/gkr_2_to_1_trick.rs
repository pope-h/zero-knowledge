use crate::{gkr_protocol::Circuit, multi_linear::MultiLinearPoly, transcript::Transcript};
use ark_ff::PrimeField;

impl<F: PrimeField> Circuit<F> {
    pub fn gkr_trick(
        &self,
        r_b: F,
        r_c: F,
        index: usize,
    ) -> (MultiLinearPoly<F>, MultiLinearPoly<F>) {
        let mut transcript = Transcript::new();

        let alpha = F::from_be_bytes_mod_order(&transcript.squeeze());
        let beta = F::from_be_bytes_mod_order(&transcript.squeeze());

        let (add_i_plus_1, mul_i_plus_1) = self.layer_i_add_mul(index + 1);

        let mut add_mle = MultiLinearPoly::new(add_i_plus_1);
        let mut mul_mle = MultiLinearPoly::new(mul_i_plus_1);

        let add_rb = add_mle.partial_evaluate(r_b, 0);
        let add_rc = add_mle.partial_evaluate(r_c, 0);

        let mul_rb = mul_mle.partial_evaluate(r_b, 0);
        let mul_rc = mul_mle.partial_evaluate(r_c, 0);

        let new_add = MultiLinearPoly::new(
            add_rb
                .computation
                .iter()
                .zip(add_rc.computation.iter())
                .map(|(&rb_val, &rc_val)| (alpha * rb_val) + (beta * rc_val))
                .collect(),
        );
        let new_mul = MultiLinearPoly::new(
            mul_rb
                .computation
                .iter()
                .zip(mul_rc.computation.iter())
                .map(|(&rb_val, &rc_val)| (alpha * rb_val) + (beta * rc_val))
                .collect(),
        );

        (new_add, new_mul)
    }
}

#[cfg(test)]
mod test {
    // use super::*;
    use ark_bn254::Fq;
    use crate::gkr_protocol::{Circuit, Gate, GateOp, Layer};

    #[test]
    fn setup_test_circuit() {
        let inputs = vec![
            Fq::from(1),
            Fq::from(2),
            Fq::from(3),
            Fq::from(4),
            Fq::from(5),
            Fq::from(6),
            Fq::from(7),
            Fq::from(8),
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

        let r_b = Fq::from(2);
        let r_c = Fq::from(3);

        circuit.add_layer(layer_1);
        circuit.add_layer(layer_2);
        circuit.add_layer(layer_3);

        let (new_add, new_mul) = circuit.gkr_trick(r_b, r_c, 0);

        dbg!(new_add);
        dbg!(new_mul);
    }
}