use sha2::{Digest, Sha256};

//=========================================================================================
// For Input: [1, 2, 3, 4, 5, 6, 7, 8]
// layers: vec![
// Layer 0 (Leaves)
// vec![H(1), H(2), H(3), H(4), H(5), H(6), H(7), H(8)],
// Layer 1
// vec![H(H(1) || H(2)), H(H(3) || H(4)), H(H(5) || H(6)), H(H(7) || H(8))],
// Layer 2
// vec![H(H(H(1) || H(2)) || H(H(3) || H(4))), H(H(H(5) || H(6)) || H(H(7) || H(8)))],
// Layer 3 (Root)
// vec![H(H(H(H(1) || H(2)) || H(H(3) || H(4))) || H(H(H(5) || H(6)) || H(H(7) || H(8))))],
// ],

// layers is a Vec<Vec<Vec<u8>>> where 
    // the outer Vec represents the layers of the Merkle tree,
    // the middle Vec represents the nodes in a layer, and 
    // the inner Vec represents the hash of a node. i.e. hash fn returns Vec<u8>
//=========================================================================================
#[derive(Debug)]
pub struct MerkleTree {
    pub layers: Vec<Vec<Vec<u8>>>,
}

//=========================================================================================
// siblings refers to the sibling nodes that are required to reconstruct the path
// from a specific leaf (input) to the root of the tree
//=========================================================================================
#[derive(Debug, Clone)]
pub struct MerkleProof {
    pub siblings: Vec<Vec<u8>>,
}

impl MerkleTree {
    pub fn new(data: &[&[u8]]) -> Self {
        let leaves = data
            .iter()
            .map(|x| MerkleTree::hash(x))
            .collect::<Vec<Vec<u8>>>();

        if leaves.is_empty() {
            return MerkleTree { layers: vec![] };
        }

        let mut layers = vec![leaves];
        while layers.last().unwrap().len() > 1 {
            let current_layer = layers.last().unwrap();
            let mut next_layer = Vec::new();

            let mut i = 0;
            while i < current_layer.len() {
                if i + 1 < current_layer.len() {
                    let combined =
                        [current_layer[i].as_slice(), current_layer[i + 1].as_slice()].concat();
                    next_layer.push(MerkleTree::hash(&combined));
                    i += 2;
                } else {
                    let combined =
                        [current_layer[i].as_slice(), current_layer[i].as_slice()].concat();
                    next_layer.push(MerkleTree::hash(&combined));
                    i += 1;
                }
            }
            layers.push(next_layer);
        }

        MerkleTree { layers }
    }

    pub fn root(&self) -> Option<Vec<u8>> {
        self.layers.last().and_then(|layer| layer.first().cloned())
    }

    pub fn generate_proof(&self, leaf: &[u8]) -> Option<MerkleProof> {
        let leaf_hash = MerkleTree::hash(leaf);
        let index = self.layers[0].iter().position(|x| x == &leaf_hash)?;

        if self.layers.is_empty() || index >= self.layers[0].len() {
            return None;
        }

        let mut siblings = Vec::new();
        let mut current_index = index;

        for layer in self.layers.iter().take(self.layers.len() - 1) {
            let sibling_index = if current_index % 2 == 0 {
                current_index + 1
            } else {
                current_index - 1
            };

            let sibling_hash = if sibling_index >= layer.len() {
                layer[current_index].clone()
            } else {
                layer[sibling_index].clone()
            };

            siblings.push(sibling_hash);
            current_index /= 2;
        }

        Some(MerkleProof { siblings })
    }

    pub fn hash(data: &[u8]) -> Vec<u8> {
        let mut hasher = Sha256::new();
        hasher.update(data);
        hasher.finalize().to_vec()
    }

    pub fn verify_proof(&self, leaf_data: &[u8], proof: &MerkleProof) -> bool {
        let leaf_hash = MerkleTree::hash(leaf_data);

        //=========================================================================================
        // Find the index of the leaf
        //=========================================================================================
        let current_index = match self.layers[0].iter().position(|x| x == &leaf_hash) {
            Some(index) => index,
            None => return false,
        };

        let root = match self.root() {
            Some(r) => r,
            None => return false,
        };

        //=========================================================================================
        // Start verification with the leaf hash
        //=========================================================================================
        let mut current_hash = leaf_hash;
        let mut idx = current_index;

        for sibling_hash in &proof.siblings {
            let combined = if idx % 2 == 0 {
                [current_hash.as_slice(), sibling_hash.as_slice()].concat()
            } else {
                [sibling_hash.as_slice(), current_hash.as_slice()].concat()
            };

            current_hash = MerkleTree::hash(&combined);
            idx /= 2;
        }

        current_hash == root
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merkle_tree() {
        let data: Vec<&[u8]> = vec![
            b"hello", b"world", b"foo", b"bar", b"baz", b"qux", b"quux", b"corge",
        ];

        let tree = MerkleTree::new(&data);

        for leaf in data.iter() {
            let proof = tree.generate_proof(leaf).unwrap();
            assert!(tree.verify_proof(leaf, &proof));
        }
    }

    #[test]
    fn test_verify_proof() {
        // Original inputs
        let data: Vec<&[u8]> = vec![
            b"hello", b"world", b"foo", b"bar", b"baz", b"qux", b"quux", b"corge",
        ];

        // Build the Merkle tree
        let tree = MerkleTree::new(&data);
        let input_to_prove = b"foo";

        // Generate the proof
        let proof = tree.generate_proof(input_to_prove).unwrap();
        println!(
            "Proof for 'foo': {:?}",
            proof
                .siblings
                .iter()
                .map(|s| hex::encode(s))
                .collect::<Vec<_>>()
        );

        // Verify the proof
        let is_valid = tree.verify_proof(input_to_prove, &proof);
        assert!(is_valid);
    }

    #[test]
    fn test_verify_proof_num() {
        let input = vec![1, 2, 3, 4, 5, 6, 7, 8];

        // Store the String values in a vector to ensure they outlive `data`
        let string_data: Vec<String> = input.iter().map(|num| num.to_string()).collect();

        // Create the `data` vector from references to the String values
        let data: Vec<&[u8]> = string_data.iter().map(|s| s.as_bytes()).collect();

        // Build the Merkle tree
        let tree = MerkleTree::new(&data);
        let input_to_prove = b"6";

        // Generate the proof
        let proof = tree.generate_proof(input_to_prove).unwrap();
        println!(
            "Proof for '6': {:?}",
            proof
                .siblings
                .iter()
                .map(|s| hex::encode(s))
                .collect::<Vec<_>>()
        );

        // Verify the proof
        let is_valid = tree.verify_proof(input_to_prove, &proof);
        assert!(is_valid);
    }
}
