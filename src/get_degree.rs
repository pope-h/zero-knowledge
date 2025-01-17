pub fn degree(dense_array: Vec<u32>) -> u32 {
    // get the length of the dense array
    (dense_array.len() - 1).try_into().unwrap()
}