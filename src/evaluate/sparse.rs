pub fn sparse_repr(x: u32, sparse_array: Vec<(u32, u32)>) -> u32 {
    let n = sparse_array.len();
    let mut result = 0;

    for i in 0..n {
        let (x_i, y_i) = sparse_array[i];
        result += x_i * x.pow(y_i);
    }

    result
}