use polynomials::*;
use shamir_secret_sharing::*;
use ark_bn254::Fq;

fn main() {
    let secret = Fq::from(4);
    let threshold = 3;
    let num_shares = 7;

    // Split the secret into shares
    let shares = generate_shares(secret, threshold, num_shares);
    println!("Shares: {:?}", &shares);

    // Can reconstruct with any threshold or more shares
    let some_shares = &shares[2..7]; // Using shares 2, 3, 4, 5, and 6
    let reconstructed = reconstruct_secret(some_shares, threshold);

    println!("Original secret: {}", secret);
    println!("Reconstructed secret: {}", reconstructed);
}
