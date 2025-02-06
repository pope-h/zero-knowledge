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

    pub fn evaluate(&self) -> Vec<F> {
        let mut current_layer = self.inputs.clone();

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

            current_layer = next_layer;
        }
        current_layer
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
        assert_eq!(result, vec![Fq::from(3), Fq::from(2)]);
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

        let result = circuit.evaluate();
        assert_eq!(result, vec![Fq::from(15)]);
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

        let result = circuit.evaluate();
        assert_eq!(result, vec![Fq::from(15), Fq::from(1680)]);
    }
}