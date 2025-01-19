use rand::Rng;

// generate any amount of points
pub fn generate_points(
    num_points: u8,
    x_max: u32,
    y_max: u32,
) -> (Vec<f64>, Vec<f64>) {
    let mut x_values: Vec<f64> = Vec::new();
    let mut y_values: Vec<f64> = Vec::new();

    while x_values.len() < num_points as usize {
        let x = rand::thread_rng().gen_range(1, x_max) as f64;
        if !x_values.contains(&x) {
            x_values.push(x);
        }
    }

    while y_values.len() < num_points as usize {
        let y = rand::thread_rng().gen_range(1, y_max) as f64;
        if !y_values.contains(&y) {
            y_values.push(y);
        }
    }

    let xs = sort_values(x_values);
    let ys = sort_values(y_values);

    // points
    (xs, ys)
}

// sort array of f64 in ascending order
pub fn sort_values(mut values: Vec<f64>) -> Vec<f64> {
    // let mut values = values;
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    values
}

// generate random x values
pub fn generate_x_values(num_points: u8, x_max: f64) -> Vec<f64> {
    let mut x_values = Vec::new();
    for _ in 0..num_points {
        let x = rand::thread_rng().gen_range(1.0, x_max);
        if !x_values.contains(&x) {
            x_values.push(x);
        }
    }
    let x_sorted = sort_values(x_values);
    x_sorted
}