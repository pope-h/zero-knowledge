use ark_ff::FftField;

// implement the FFT algorithm
pub struct FastFourierTransform<F: FftField> {
    pub coefficients: Vec<F>,
}

impl<F: FftField> FastFourierTransform<F> {
    fn new(coefficients: Vec<F>) -> Self {
        FastFourierTransform { coefficients }
    }

    fn is_power_of_two(&self) -> bool {
        self.coefficients.len().is_power_of_two()
    }

    // this splits the array by the even and odd indexes
    fn split_array(&self) -> (Vec<F>, Vec<F>) {
        let mut even = Vec::new();
        let mut odd = Vec::new();

        for (i, coeff) in self.coefficients.iter().enumerate() {
            if i % 2 == 0 {
                even.push(*coeff);
            } else {
                odd.push(*coeff);
            }
        }

        (even, odd)
    }

    fn root_of_unity(&self) -> F {
        F::get_root_of_unity(self.coefficients.len() as u64).unwrap()
    }

    // this is the FFT function i.e. converting Coeff => Values
    pub fn evaluate(&self) -> Self {
        if !self.is_power_of_two() {
            panic!("The computation array must be in the power of 2");
        }

        let n = self.coefficients.len() as u64;
        let w = self.root_of_unity();

        if n == 1 {
            return FastFourierTransform {
                coefficients: self.coefficients.clone(),
            };
        }

        let (even, odd) = self.split_array();

        let y_even = FastFourierTransform::new(even).evaluate();
        let y_odd = FastFourierTransform::new(odd).evaluate();

        let mut y = vec![F::zero(); n as usize];

        for i in 0..(n / 2) {
            y[i as usize] = y_even.coefficients[i as usize]
                + (w.pow(&[i]) * y_odd.coefficients[i as usize]);
            y[(i + n / 2) as usize] = y_even.coefficients[i as usize]
                - (w.pow(&[i]) * y_odd.coefficients[i as usize]);
        }

        FastFourierTransform { coefficients: y }
    }

    // this is the inverse IFFT function i.e. converting Values => Coeff
    pub fn interpolate(&self) -> Self {
        if !self.is_power_of_two() {
            panic!("The computation array must be in the power of 2");
        }

        let n = self.coefficients.len() as u64;
        let w = self.root_of_unity().inverse().unwrap();
        if n == 1 {
            return FastFourierTransform {
                coefficients: self.coefficients.clone(),
            };
        }

        let (even, odd) = self.split_array();

        let y_even = FastFourierTransform::new(even).interpolate();
        let y_odd = FastFourierTransform::new(odd).interpolate();

        let mut y = vec![F::zero(); n as usize];

        for i in 0..(n / 2) {
            y[i as usize] =
                y_even.coefficients[i as usize] + (w.pow(&[i]) * y_odd.coefficients[i as usize]);
            y[(i + n / 2) as usize] =
                y_even.coefficients[i as usize] - (w.pow(&[i]) * y_odd.coefficients[i as usize]);
        }
        
        // let y_divided: Vec<F> = y.iter_mut().map(|elem| *elem / F::from(self.coefficients.len() as u64)).collect();

        FastFourierTransform { coefficients: y }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_377::Fr;

    #[test]
    fn test_fft() {
        let coefficients = vec![Fr::from(5), Fr::from(0), Fr::from(0), Fr::from(2)];
        let undivided_coefficients = vec![Fr::from(20), Fr::from(0), Fr::from(0), Fr::from(8)];

        let fft = FastFourierTransform::new(coefficients.clone());
        let values = fft.evaluate();
        let interpolated = values.interpolate();

        // assert_eq!(interpolated.coefficients, coefficients);
        assert_eq!(interpolated.coefficients, undivided_coefficients);
    }
}
