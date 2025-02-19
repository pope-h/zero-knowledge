use ark_ff::{BigInteger, PrimeField};

#[derive(Debug, PartialEq, Clone)]
pub struct MultiLinearPoly<F: PrimeField> {
    // 2a + 3b
    // computation is simply passing in 00, 01, 10, 11 and getting the result in an array
    // 00 -> 0, 01 -> 3, 10 -> 2, 11 -> 5 => [0, 3, 2, 5]
    pub computation: Vec<F>,
}

impl<F: PrimeField> MultiLinearPoly<F> {
    pub fn new(computation: Vec<F>) -> Self {
        if !computation.len().is_power_of_two() {
            panic!("The computation array must be in the power of 2");
        }
        MultiLinearPoly { computation }
    }

    fn variable_count(&self) -> u32 {
        // if the variable is e.g. a, b = 2
        // the computation array will have a size of 2^2 = 4
        // thus the variable count is log2 of the length of the computation array
        self.computation.len().ilog2()
    }

    fn get_power(&self, eval_point_index: usize) -> usize {
        self.variable_count() as usize - eval_point_index - 1
    }

    pub fn partial_evaluate(&mut self, eval_value: F, eval_value_position: usize) -> Self {
        let new_length = self.computation.len() / 2;
        let mut new_computation = Vec::with_capacity(new_length);

        let step = 2u32.pow(self.get_power(eval_value_position) as u32);

        let mut y_1_i = 0;

        while y_1_i < self.computation.len() {
            let y_1 = self.computation[y_1_i];
            let y_2 = self.computation[y_1_i + step as usize];

            new_computation.push(y_1 + (y_2 - y_1) * eval_value);

            if (y_1_i + 1) % step as usize == 0 {
                y_1_i = y_1_i + step as usize + 1;
            } else {
                y_1_i += 1;
            }
        }

        MultiLinearPoly::new(new_computation)
    }

    pub fn evaluate(&mut self, eval_points: Vec<F>) -> Self {
        if eval_points.len() != self.variable_count() as usize {
            panic!("The number of eval points must be equal to the number of variables");
        }

        let mut this_computation = MultiLinearPoly::new(self.computation.clone());
        let mut i = 0;

        while i < eval_points.len() {
            this_computation = this_computation.partial_evaluate(eval_points[i], 0);
            i += 1;
        }

        this_computation
    }

    pub fn to_bytes(computation: Vec<F>) -> Vec<u8> {
        computation
            .iter()
            .flat_map(|x| F::into_bigint(*x).to_bytes_be())
            .collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    pub fn setup_mle_poly() -> MultiLinearPoly<Fq> {
        let computation = vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(4),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(3),
            Fq::from(5),
            Fq::from(9),
            Fq::from(8),
            Fq::from(12),
        ];
        let multi_linear_poly = MultiLinearPoly::new(computation);

        multi_linear_poly
    }

    #[test]
    fn test_partial_evaluate() {
        // 2a + 3b
        let computation = vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)];
        let mut multi_linear_poly = MultiLinearPoly::new(computation);

        let eval_value = Fq::from(1);
        let eval_value_position = 0;
        let result = multi_linear_poly.partial_evaluate(eval_value, eval_value_position);

        assert_eq!(result.computation, vec![Fq::from(2), Fq::from(5)]);

        let eval_value_b = Fq::from(4);
        let eval_value_b_position = 1;
        let result_b = multi_linear_poly.partial_evaluate(eval_value_b, eval_value_b_position);

        assert_eq!(result_b.computation, vec![Fq::from(12), Fq::from(14)]);
    }

    #[test]
    fn test_partial_evaluate_at_last_index() {
        // 00030025 0 9 0 11 at c = 3
        let computation = vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(3),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(5),
        ];
        let mut multi_linear_poly = MultiLinearPoly::new(computation);

        let eval_value = Fq::from(3);
        let eval_value_position = 2;
        let result = multi_linear_poly.partial_evaluate(eval_value, eval_value_position);

        assert_eq!(
            result.computation,
            vec![Fq::from(0), Fq::from(9), Fq::from(0), Fq::from(11)]
        );
    }

    #[test]
    fn test_evaluate() {
        let computation = vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)];
        let mut multi_linear_poly = MultiLinearPoly::new(computation);

        let eval_points = vec![Fq::from(1), Fq::from(1)];
        let result = multi_linear_poly.evaluate(eval_points);

        assert_eq!(result.computation, vec![Fq::from(5)]);
    }

    #[test]
    fn test_partial_evaluate_2() {
        let mut multi_linear_poly = setup_mle_poly();

        let eval_value = Fq::from(4);
        let eval_value_position = 0;
        let result = multi_linear_poly.partial_evaluate(eval_value, eval_value_position);

        assert_eq!(
            result.computation,
            vec![
                Fq::from(0),
                Fq::from(0),
                Fq::from(12),
                Fq::from(12),
                Fq::from(20),
                Fq::from(24),
                Fq::from(32),
                Fq::from(36)
            ]
        );
    }

    #[test]
    fn test_partial_evaluate_2_at_2_points() {
        let mut multi_linear_poly = setup_mle_poly();
        let eval_value_position = 0;

        let eval_value_a = Fq::from(4);

        let eval_value_b = Fq::from(2);
        let result = multi_linear_poly
            .partial_evaluate(eval_value_a, eval_value_position)
            .partial_evaluate(eval_value_b, eval_value_position);

        assert_eq!(
            result.computation,
            vec![Fq::from(40), Fq::from(48), Fq::from(52), Fq::from(60)]
        );
    }

    #[test]
    fn test_evaluate_2() {
        let mut multi_linear_poly = setup_mle_poly();

        let eval_points = vec![Fq::from(4), Fq::from(2), Fq::from(6), Fq::from(1)];
        let result = multi_linear_poly.evaluate(eval_points);

        assert_eq!(result.computation, vec![Fq::from(120)]);
    }

    #[test]
    fn test_to_bytes() {
        let computation = vec![Fq::from(5)];
        let bytes = MultiLinearPoly::to_bytes(computation);

        assert_eq!(
            bytes,
            vec![
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 5
            ]
        );
    }

    // Test for understanding GKR
    #[test]
    fn test_mul_add_evaluate() {
        let computation = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)];
        let mut poly = MultiLinearPoly::new(computation);

        // using a random challenge r = 2 for a
        let eval_value = Fq::from(0);
        let eval_position = 0;
        let result = poly.partial_evaluate(eval_value, eval_position);

        dbg!(result);
        // 21888242871839275222246405745257275088696311157297823662689037894645226208582
    }

    #[test]
    fn test_w_evaluate() {
        let computation = vec![Fq::from(3), Fq::from(7), Fq::from(11), Fq::from(56)];
        let mut poly = MultiLinearPoly::new(computation);

        let eval_points = vec![Fq::from(1), Fq::from(1)];
        let result = poly.evaluate(eval_points);

        assert_eq!(result.computation, vec![Fq::from(56)]);
    }

    #[test]
    fn test_w_evaluate2() {
        let computation = vec![Fq::from(0), Fq::from(12)];
        let mut poly = MultiLinearPoly::new(computation);

        let eval_points = vec![Fq::from(2)];
        let result = poly.evaluate(eval_points);

        dbg!(result);
    }
}
