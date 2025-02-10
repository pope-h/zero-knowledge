use ark_ff::PrimeField;

#[derive(Debug, PartialEq)]
pub enum GateOp {
    Add,
    Mul,
}

pub struct Gate {
    pub left: usize,
    pub right: usize,
    pub op: GateOp,
    pub output: usize,
}

pub struct Layer {
    pub gates: Vec<Gate>,
}

pub struct Circuit<F: PrimeField> {
    pub inputs: Vec<F>,
    pub layers: Vec<Layer>,
}

impl<F: PrimeField> Circuit<F> {
    pub fn new(inputs: Vec<F>) -> Self {
        Circuit { layers: Vec::new(), inputs }
    }

    pub fn add_layer(&mut self, layer: Layer) {
        self.layers.push(layer);
    }

    pub fn evaluate(&self) -> Vec<Vec<F>> {
        let mut current_layer = self.inputs.clone();
        let mut eval_layers = vec![];

        for layer in self.layers.iter() {
            let mut next_layer = vec![F::zero(); layer.gates.len()];

            for gate in layer.gates.iter() {
                let left = current_layer[gate.left];
                let right = current_layer[gate.right];

                let result = match gate.op {
                    GateOp::Add => left + right,
                    GateOp::Mul => left * right,
                };

                next_layer[gate.output] = result;
            }

            eval_layers.push(next_layer.clone());
            current_layer = next_layer;
        }
        eval_layers
    }

    pub fn compute_layer_i_add_mul(&self, layer_i: usize) -> (Vec<F>, Vec<F>) {
        let layer = &self.layers[self.layers.len() - layer_i - 1];

        // my discovery is that for a 2-gate with 1-bit each at output and 2-bits each at input
        // for gate counts not in the power of 2, we pad it to the next power of 2
        let gate_count;
        let input_bits;
        let output_bits;

        if layer.gates.len() == 1 {
            input_bits = 1;
            output_bits = 1;
        } else if layer.gates.len().is_power_of_two() {
            gate_count = layer.gates.len();

            input_bits = (2 * gate_count).ilog2();
            output_bits = gate_count.ilog2();
        } else {
            gate_count = layer.gates.len().next_power_of_two();

            input_bits = (2 * gate_count).ilog2();
            output_bits = gate_count.ilog2();
        }

        let n_bits = (2 * input_bits) + output_bits;    // (left_bit + right_bit )+ output_bit
        let total_combinations = 2usize.pow(n_bits);
        let mut add_vec = vec![F::zero(); total_combinations];
        let mut mul_vec = vec![F::zero(); total_combinations];

        let output_start_index = n_bits - output_bits;
        let left_start_index = output_start_index - input_bits;
        let right_start_index = 0 as u32;

        for gate in &layer.gates {
            let index = match gate.op {
                GateOp::Mul => (gate.output << output_start_index) + (gate.left << left_start_index) + (gate.right << right_start_index),
                GateOp::Add => (gate.output << output_start_index) + (gate.left << left_start_index) + (gate.right << right_start_index),
            };

            match gate.op {
                GateOp::Mul => mul_vec[index] = F::one(),
                GateOp::Add => add_vec[index] = F::one(),
            }
        }

        (add_vec, mul_vec)
    }
    
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    fn test_gate_creation() {
        let gate = Gate { left: 0, right: 1, op: GateOp::Add, output: 0 };
        assert_eq!(gate.left, 0);
        assert_eq!(gate.right, 1);
        assert_eq!(gate.op, GateOp::Add);
        assert_eq!(gate.output, 0);
    }

    #[test]
    fn test_layer_creation() {
        let gate_1 = Gate { left: 0, right: 1, op: GateOp::Add, output: 0 };
        let gate_2 = Gate { left: 0, right: 1, op: GateOp::Mul, output: 1 };
        let layer = Layer { gates: vec![gate_1, gate_2] };

        assert_eq!(layer.gates.len(), 2);
        assert_eq!(layer.gates[0].op, GateOp::Add);
        assert_eq!(layer.gates[1].op, GateOp::Mul);
    }

    #[test]
    fn test_evaluate() {
        let inputs = vec![Fq::from(1), Fq::from(2)];
        let mut circuit = Circuit::new(inputs);

        let layer = Layer {
            gates: vec![
                Gate { left: 0, right: 1, op: GateOp::Add, output: 0 },
                Gate { left: 0, right: 1, op: GateOp::Mul, output: 1 },
            ]
        };

        circuit.add_layer(layer);

        let result = circuit.evaluate();
        assert_eq!(result, vec![vec![Fq::from(3), Fq::from(2)]]);
    }

    #[test]
    fn test_evaluate_2() {
        let inputs = vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)];
        let mut circuit = Circuit::new(inputs);

        let layer_1 = Layer {
            gates: vec![
                Gate { left: 0, right: 1, op: GateOp::Add, output: 0 },
                Gate { left: 2, right: 3, op: GateOp::Mul, output: 1 },
            ]
        };

        let layer_2 = Layer {
            gates: vec![
                Gate { left: 0, right: 1, op: GateOp::Add, output: 0 },
            ]
        };

        circuit.add_layer(layer_1);
        circuit.add_layer(layer_2);

        let eval_layer_1 = vec![Fq::from(3), Fq::from(12)];
        let eval_layer_2 = vec![Fq::from(15)];

        let result = circuit.evaluate();
        assert_eq!(result, vec![eval_layer_1, eval_layer_2]);
    }

    #[test]
    fn test_two_evaluate_returned() {
        let inputs = vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4), Fq::from(5), Fq::from(6), Fq::from(7), Fq::from(8)];
        let mut circuit = Circuit::new(inputs);

        let layer_1 = Layer {
            gates: vec![
                Gate { left: 0, right: 1, op: GateOp::Add, output: 0 },
                Gate { left: 2, right: 3, op: GateOp::Mul, output: 1 },
                Gate { left: 4, right: 5, op: GateOp::Mul, output: 2 },
                Gate { left: 6, right: 7, op: GateOp::Mul, output: 3 },
            ]
        };

        let layer_2 = Layer {
            gates: vec![
                Gate { left: 0, right: 1, op: GateOp::Add, output: 0 },
                Gate { left: 2, right: 3, op: GateOp::Mul, output: 1 },
            ]
        };

        circuit.add_layer(layer_1);
        circuit.add_layer(layer_2);

        let eval_layer_1 = vec![Fq::from(3), Fq::from(12), Fq::from(30), Fq::from(56)];
        let eval_layer_2 = vec![Fq::from(15), Fq::from(1680)];

        let result = circuit.evaluate();
        assert_eq!(result, vec![eval_layer_1, eval_layer_2]);
    }

    #[test]
    fn test_arithmetic_circuit() {
        let inputs = vec![Fq::from(1), Fq::from(2), Fq::from(1), Fq::from(4)];
        let mut circuit = Circuit::new(inputs);

        let layer_1 = Layer {
            gates: vec![
                Gate { left: 0, right: 0, op: GateOp::Mul, output: 0 },
                Gate { left: 1, right: 1, op: GateOp::Mul, output: 1 },
                Gate { left: 1, right: 2, op: GateOp::Mul, output: 2 },
                Gate { left: 3, right: 3, op: GateOp::Mul, output: 3 },
            ]
        };

        let layer_2 = Layer {
            gates: vec![
                Gate { left: 0, right: 1, op: GateOp::Mul, output: 0 },
                Gate { left: 2, right: 3, op: GateOp::Mul, output: 1 },
            ]
        };

        circuit.add_layer(layer_1);
        circuit.add_layer(layer_2);

        let eval_layer_1 = vec![Fq::from(1), Fq::from(4), Fq::from(2), Fq::from(16)];
        let eval_layer_2 = vec![Fq::from(4), Fq::from(32)];

        let result = circuit.evaluate();
        assert_eq!(result, vec![eval_layer_1, eval_layer_2]);
    }

    #[test]
    fn test_compute_mul_add() {
        let inputs = vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4), Fq::from(5), Fq::from(6), Fq::from(7), Fq::from(8)];
        let mut circuit = Circuit::new(inputs);

        let layer_1 = Layer {
            gates: vec![
                Gate { left: 0, right: 1, op: GateOp::Add, output: 0 },
                Gate { left: 2, right: 3, op: GateOp::Mul, output: 1 },
                Gate { left: 4, right: 5, op: GateOp::Mul, output: 2 },
                Gate { left: 6, right: 7, op: GateOp::Mul, output: 3 },
            ]
        };

        let layer_2 = Layer {
            gates: vec![
                Gate { left: 0, right: 1, op: GateOp::Add, output: 0 },
                Gate { left: 2, right: 3, op: GateOp::Mul, output: 1 },
            ]
        };

        let layer_3 = Layer {
            gates: vec![
                Gate { left: 0, right: 1, op: GateOp::Add, output:0 }
            ]
        };

        circuit.add_layer(layer_1);
        circuit.add_layer(layer_2);
        circuit.add_layer(layer_3);

        let output = circuit.compute_layer_i_add_mul(0);
        dbg!(output);
        // assert_eq!(output, ([Fq::from(0), Fq::from(1), Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(0)], [Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(0)]));
    }
}