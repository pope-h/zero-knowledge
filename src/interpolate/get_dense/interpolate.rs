use crate::interpolate::expansion;

pub fn interpolate_dense(points: Vec<(f64, f64)>) -> Vec<f64> {
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
            numerator = expansion::multiply_poly(&numerator, &term);
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