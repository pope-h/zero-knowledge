use crate::{gkr_protocol::Circuit, multi_linear::MultiLinearPoly, transcript::Transcript};
use ark_ff::PrimeField;

impl<F: PrimeField> Circuit<F> {
    pub fn gkr_trick(
        &self,
        challenges: Vec<F>,
        index: usize,
    ) -> (MultiLinearPoly<F>, MultiLinearPoly<F>) {
        let mut transcript = Transcript::new();

        let alpha = F::from_be_bytes_mod_order(&transcript.squeeze());
        let beta = F::from_be_bytes_mod_order(&transcript.squeeze());

        let (add_i, mul_i) = self.layer_i_add_mul(index);

        let mut add_rb = MultiLinearPoly::new(add_i.clone());
        let mut add_rc = MultiLinearPoly::new(add_i);
        let mut mul_rb = MultiLinearPoly::new(mul_i.clone());
        let mut mul_rc = MultiLinearPoly::new(mul_i);

        let mid = challenges.len() / 2;
        let (r_b_challenges, r_c_challenges) = challenges.split_at(mid);

        for &r_b in r_b_challenges {
            add_rb = add_rb.partial_evaluate(r_b, 0);
            mul_rb = mul_rb.partial_evaluate(r_b, 0);
        }

        for &r_c in r_c_challenges {
            add_rc = add_rc.partial_evaluate(r_c, 0);
            mul_rc = mul_rc.partial_evaluate(r_c, 0);
        }

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

    pub fn new_claimed_sum(&self, w_i_arr: Vec<F>, challenges: Vec<F>) -> F {
        let mut transcript = Transcript::new();

        let w_i_eval = MultiLinearPoly::new(w_i_arr);

        let alpha = F::from_be_bytes_mod_order(&transcript.squeeze());
        let beta = F::from_be_bytes_mod_order(&transcript.squeeze());

        let mut w_i_b = w_i_eval.clone();
        let mut w_i_c = w_i_eval;

        let mid = challenges.len() / 2;
        let (r_b_challenges, r_c_challenges) = challenges.split_at(mid);

        for &r_b in r_b_challenges {
            w_i_b = w_i_b.partial_evaluate(r_b, 0);
        }

        for &r_c in r_c_challenges {
            w_i_c = w_i_c.partial_evaluate(r_c, 0);
        }

        // claimed_sum = (alpha * w_i(*b)) + (beta * w_i(*c))
        let claimed_sum = w_i_b
            .computation
            .iter()
            .zip(w_i_c.computation.iter())
            .map(|(&w_i_b_val, &w_i_c_val)| (alpha * w_i_b_val) + (beta * w_i_c_val))
            .sum();

        claimed_sum
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
        let challenges = vec![r_b, r_c];

        circuit.add_layer(layer_1);
        circuit.add_layer(layer_2);
        circuit.add_layer(layer_3);

        let (new_add, new_mul) = circuit.gkr_trick(challenges, 0);

        dbg!(new_add);
        dbg!(new_mul);
    }

    #[test]
    fn test_new_claimed_sum() {
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
        let challenges = vec![r_b, r_c];

        circuit.add_layer(layer_1);
        circuit.add_layer(layer_2);
        circuit.add_layer(layer_3);

        let evaluated_circuit = circuit.evaluate();
        let w_i_eval = evaluated_circuit[1].clone();

        let claimed_sum = circuit.new_claimed_sum(w_i_eval, challenges);

        dbg!(claimed_sum);
    }
}