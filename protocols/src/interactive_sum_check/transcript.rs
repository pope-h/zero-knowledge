use sha3::{Digest, Keccak256};

pub struct Transcript {
    hasher: Keccak256, // Keep the hasher as part of the state
}

impl Transcript {
    pub fn new() -> Self {
        Transcript {
            hasher: Keccak256::new(),
        }
    }

    pub fn append(&mut self, byte_array: &[u8]) {
        self.hasher.update(byte_array);
    }

    pub fn challenge(&mut self) -> Vec<u8> {
        // Clone the current hasher state to not affect future updates
        self.hasher.clone().finalize().to_vec()
    }
}
