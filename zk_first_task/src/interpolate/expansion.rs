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