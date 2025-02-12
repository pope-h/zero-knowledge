use ark_ff::PrimeField;
use crate::multi_linear::MultiLinearPoly;

/*
    * ProductPoly is a struct that represents a product of two polynomials.
    * It has a field poly_array that is a vector of MultiLinearPoly.
    * It has a method new that takes a vector of MultiLinearPoly and returns a ProductPoly.
    * It has a method get_degree that returns the length of the poly_array.
    * It has a method partial_evaluate that takes an eval_value and eval_value_position and returns a MultiLinearPoly.
    * It has a method evaluate that returns a vector of F.

    * This works such that if two polys [0, 0, 0, 2] and [0, 0, 0, 3] are passed representing 2ab and 3ab respectively,
    * Partial evaluation at 0, 1, and 2 will return [[0, 0], [0, 0]], [[0, 2], [0, 3]] and [[0, 4], [0, 6]] respectively.
       * Partial evaluate will carry out a product for each array for example [[0, 4], [0, 6]] => [0*0, 4*6] => [0, 24]
    * Evaluate gets [0, 0], [0, 6], and [0, 24] for each call and returns [0, 6, 24] after performing a sum on each array.
 */

pub struct ProductPoly<F: PrimeField> {
    pub poly_array: Vec<MultiLinearPoly<F>>,
}

impl<F: PrimeField> ProductPoly<F> {
    pub fn new(poly_array: Vec<MultiLinearPoly<F>>) -> Self {
        if poly_array[0].computation.len() != poly_array[1].computation.len() {
            panic!("The polynomials must have the same length");
        }
        ProductPoly { poly_array }
    }

    fn get_degree(&self) -> usize {
        self.poly_array.len()
    }

    pub fn partial_evaluate(&mut self, eval_value: F, eval_value_position: usize) -> MultiLinearPoly<F> {
        let size = self.poly_array[0].computation.len() / 2;
        let mut new_array: Vec<F> = Vec::with_capacity(size);
        let new_poly_array: Vec<MultiLinearPoly<F>> = self.poly_array.iter_mut().map(|m_poly| m_poly.partial_evaluate(eval_value, eval_value_position)).collect();
        // e.g. [[1, 2], [3, 4], [5, 6]] -> [[1 * 3 * 5], [2 * 4 * 6]] -> [15, 48]
        for i in 0..new_poly_array[0].computation.len() {
            let prod = new_poly_array.iter().map(|array| {
                array.computation[i]
            }).product();

            new_array.push(prod);
        }

        MultiLinearPoly::new(new_array)
    }

    pub fn evaluate(&mut self) -> Vec<F> {
        // 1 is added to the degree to satisfy the (d+1) return for points
        // for example if the degree is 2, the points will be [0, 1, 2] or if the degree is 3, the points will be [0, 1, 2, 3]
        let count = self.get_degree() + 1;
        let mut new_array = Vec::with_capacity(count);

        let mut this_computation = ProductPoly::new(self.poly_array.clone());

        let mut i = 0;
        while i < count {
            let eval_point = F::from(i as u64);
            let partial_eval = this_computation.partial_evaluate(eval_point, 0);
            let element_sum: F = partial_eval.computation.iter().sum();
            new_array.push(element_sum);
            i += 1;
        }

        new_array
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    fn test_product_poly() {
        let poly_1 = MultiLinearPoly::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)]);
        let poly_2 = MultiLinearPoly::new(vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)]);
        let mut product_poly = ProductPoly::new(vec![poly_1, poly_2]);
        let result = product_poly.evaluate();
        assert_eq!(result, vec![Fq::from(0), Fq::from(6), Fq::from(24)]);
    }

    #[test]
    fn test_product_poly2() {
        let poly_1 = MultiLinearPoly::new(vec![Fq::from(0), Fq::from(8)]);
        let poly_2 = MultiLinearPoly::new(vec![Fq::from(0), Fq::from(12)]);
        let mut product_poly = ProductPoly::new(vec![poly_1, poly_2]);
        let result = product_poly.evaluate();
        assert_eq!(result, vec![Fq::from(0), Fq::from(96), Fq::from(384)]);
    }
}