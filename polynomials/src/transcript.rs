use sha3::{Digest, Keccak256};

pub struct Transcript {
    hasher: Keccak256,  // Keep the hasher as part of the state
}

impl Transcript {
    pub fn new() -> Self {
        Transcript {
            hasher: Keccak256::new()
        }
    }

    pub fn absorb(&mut self, byte_array: &[u8]) {
        self.hasher.update(byte_array);
    }

    pub fn squeeze(&mut self) -> Vec<u8> {
        // Clone the current hasher state to not affect future updates
        let hasher_clone = self.hasher.clone();
        let challenge_hash = hasher_clone.finalize().to_vec();
        self.hasher.update(challenge_hash.clone());
    
        challenge_hash
    }
}