use rand::Rng;
use polynomials::*;
use shamir_secret_sharing::*;

fn generate_threshold() -> u8 {
    rand::thread_rng().gen_range(2, 10)
}

fn main() {
    // Generate random parameters
    let num_points = generate_threshold();
    let x_max = 10;
    let y_max = 20;

    println!("Generating {} points with x_max = {}, y_max = {}", num_points, x_max, y_max);

    // Generate points
    let (xs, ys) = generate_points(num_points, x_max, y_max);
    
    println!("\nGenerated Points:");
    for (i, (x, y)) in xs.iter().zip(ys.iter()).enumerate() {
        println!("Point {}: ({:.2}, {:.2})", i + 1, x, y);
    }

    // Create polynomial
    let poly = UnivariatePoly::interpolate(xs.clone(), ys.clone());
    
    println!("\nPolynomial Details:");
    println!("Degree: {}", poly.degree());
    println!("Coefficients: {:.2?}", poly.coefficient);

    // Test evaluation at a few points
    println!("\nTesting polynomial evaluation:");
    for (x, y) in xs.iter().zip(ys.iter()) {
        let evaluated = poly.evaluate(*x);
        println!("f({:.2}) = {:.2} (original y = {:.2})", x, evaluated, y);
    }
}