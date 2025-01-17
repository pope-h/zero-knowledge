use crate::interpolate::expansion;

pub fn interpolate_sparse(points: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    let n = points.len();
    let mut coefficients = vec![0.0; n];
    
    // Calculate coefficients
    for i in 0..n {
        let (x_i, y_i) = points[i];
        let mut numerator = vec![1.0];
        let mut denominator = 1.0;
        
        for j in 0..n {
            if i == j {
                continue;
            }
            let (x_j, _) = points[j];
            
            let initial_coeffs = vec![-x_j, 1.0];
            numerator = expansion::multiply_poly(&numerator, &initial_coeffs);
            denominator *= x_i - x_j;
        }
        
        let multiplier = y_i / denominator;
        for k in 0..numerator.len() {
            if k < n {
                coefficients[k] += numerator[k] * multiplier;
            }
        }
    }
    
    let mut final_result = Vec::new();
    
    // Switch the order in the tuple to (coefficient, power)
    for (power, &coeff) in coefficients.iter().enumerate() {
        if coeff.abs() > 1e-10 {
            final_result.push((coeff, power as f64));  // Swapped order here
        }
    }
    
    // Sort in descending order by coefficient
    final_result.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    
    final_result
}