pub fn evaluate(x: u32 ,dense_array: Vec<u32>) -> u32 {
    // returns the value of the polynomial at x
    let mut new_array = Vec::new();
    for (i, val) in dense_array.iter().enumerate() {
        new_array.push(val * x.pow(i.try_into().unwrap()));
    }

    new_array.iter().sum()
    
}

pub fn degree(dense_array: Vec<u32>) -> u32 {
    // get the length of the dense array
    (dense_array.len() - 1).try_into().unwrap()
}

/*  
    This function expands any two polynomials at a time i.e. y = (x - a)(x - b) is done in following steps 
    y += x * x; result += x * -b; result += -a * x; result += -a * -b
    then result is returned
*/
pub fn multiply_poly(p1: &[f64], p2: &[f64]) -> Vec<f64> {
    let n1 = p1.len();
    let n2 = p2.len();
    let mut result = vec![0.0; n1 + n2 - 1];
    
    for i in 0..n1 {
        for j in 0..n2 {
            result[i + j] += p1[i] * p2[j];
        }
    }
    result
}

pub fn interpolate(points: Vec<(f64, f64)>) -> Vec<f64> {
    let n = points.len();
    let mut final_result = vec![0.0; n];
    
    // For each point
    for i in 0..n {
        let (x_i, y_i) = points[i];
        let mut numerator = vec![1.0]; // Start with 1
        let mut denominator = 1.0;
        
        // Build the Lagrange basis polynomial
        for j in 0..n {
            if i == j { continue; }
            let (x_j, _) = points[j];
            
            // Create (x - x_j) term
            let term = vec![-x_j, 1.0];  // coefficients of (-x_j + 1.0x) is taken. Recall that -x_j is a constant hence why it is used in full
            numerator = multiply_poly(&numerator, &term);
            denominator *= x_i - x_j;
        }
        
        // Scale each coefficient by y_i/denominator
        let scale = y_i / denominator;
        for k in 0..numerator.len() {
            if k < n {
                final_result[k] += numerator[k] * scale;
            }
        }
    }
    
    final_result
}

#[cfg(test)]
mod tests {
    // use super::*;

    #[test]
    fn test_evaluate() {
        let dense_array = vec![5, 0, 0, 2];
        let x = 2;
        let result = super::evaluate(x, dense_array);
        assert_eq!(result, 21);
    }

    #[test]
    fn test_degree() {
        let dense_array = vec![5, 0, 0, 2];
        let dense_array_degree = super::degree(dense_array);
        assert_eq!(dense_array_degree, 3);
    }

    #[test]
    fn test_interpolate() {
        let points = vec![(0.0, 5.0), (1.0, 7.0), (2.0, 21.0), (3.0, 59.0)];
        let result = super::interpolate(points);
        assert_eq!(result, vec![5.0, 0.0, 0.0, 2.0]);
    }
}
