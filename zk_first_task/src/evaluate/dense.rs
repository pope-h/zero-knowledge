pub fn dense_repr(x: u32 ,dense_array: Vec<u32>) -> u32 {
    // returns the value of the polynomial at x
    let mut new_array = Vec::new();
    for (i, val) in dense_array.iter().enumerate() {
        new_array.push(val * x.pow(i.try_into().unwrap()));
    }

    new_array.iter().sum()
    
}