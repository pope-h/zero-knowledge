pub mod get_degree;
pub mod evaluate;
pub mod interpolate;

#[cfg(test)]
mod tests {
    // use super::*;

    use crate::get_degree;
    use crate::evaluate::{dense, sparse};
    use crate::interpolate::get_dense::interpolate;

    #[test]
    fn test_dense_evaluate() {
        let dense_array = vec![5, 0, 0, 2];
        let x = 2;
        let result = dense::dense_repr(x, dense_array);
        assert_eq!(result, 21);
    }

    #[test]
    fn test_sparse_evaluate() {
        let sparse_array = vec![(2, 3), (5, 0)];
        let x = 2;
        let result = sparse::sparse_repr(x, sparse_array);
        assert_eq!(result, 21);
    }

    #[test]
    fn test_degree() {
        let dense_array = vec![5, 0, 0, 2];
        let dense_array_degree = get_degree::degree(dense_array);
        assert_eq!(dense_array_degree, 3);
    }

    #[test]
    fn test_interpolate() {
        let points = vec![(0.0, 5.0), (1.0, 7.0), (2.0, 21.0), (3.0, 59.0)];
        let result = interpolate::interpolate_dense(points);
        assert_eq!(result, vec![5.0, 0.0, 0.0, 2.0]);
    }
}
