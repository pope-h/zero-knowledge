pub mod shamir_secret_sharing;

use std::iter::{Product, Sum};
use std::ops::{Add, Mul};

#[derive(Debug, PartialEq, Clone)]
pub struct UnivariatePoly {
    // 1 coefficient for each power of x
    pub coefficient: Vec<f64>,
}

impl UnivariatePoly {
    fn new(coefficient: Vec<f64>) -> Self {
        UnivariatePoly { coefficient }
    }

    pub fn degree(&self) -> usize {
        self.coefficient.len() - 1
    }

    pub fn evaluate(&self, x: f64) -> f64 {
        // self.coefficient
        //     .iter()
        //     .enumerate()
        //     .map(|(i, coeff)| coeff * x.powf(i as f64))
        //     .sum()

        /*
            rev() is used here cos the lowest power of x is the first element in the inputted vector e.g. 2x^2 + 3x + 1 => [1, 3, 2]
            This is then reversed to [2, 3, 1] so that the highest power of x is the first element in the vector
            We can then run 2 * 2 + 3 = 7; 7 * 2 + 1 = 15 evaluating at x = 2;

            Broader explanation:
            with 5x^4 + 3x^2 + 7x + 11 = 5x^4 + 0x^3 + 3x^2 + 7x + 11
            (((5x + 0)x + 3)x + 7)x + 11

            in the program [11, 7, 3, 0, 5] will be inputted but this is then reversed to [5, 0, 3, 7, 11]
            (((5x + 0)x + 3)x + 7)x + 11 = 5*2 + 0 = 10; 10*2 + 3 = 23; 23*2 + 7 = 53; 53*2 + 11 = 117 evaluating at x = 2
        */
        self.coefficient
            .iter()
            .rev()
            .cloned()
            .reduce(|acc, curr| acc * x + curr)
            .unwrap()
    }

    pub fn interpolate(xs: Vec<f64>, ys: Vec<f64>) -> Self {
        xs.iter()
            .zip(ys.iter())
            .map(|(x, y)| Self::basis(x, &xs).scalar_mul(y))
            .sum()
    }

    // Multiplies an array with an integer
    fn scalar_mul(&self, scalar: &f64) -> Self {
        UnivariatePoly {
            coefficient: self.coefficient.iter().map(|x| x * scalar).collect(),
        }
    }

    /*
        [1, 2, 3]
        L_2(x) = (x - 1)(x - 3)
                 --------------
                 (2 - 1)(2 - 3)
        x - 1
        [-x, 1]

        [1, 2, 3] -> [1, 3] -> [(x - 1), (x - 3)]
    */
    fn basis(x: &f64, interpolating_set: &[f64]) -> Self {
        //numerator
        let numerator: UnivariatePoly = interpolating_set
            .iter()
            .filter(|val| *val != x)
            .map(|x_i| UnivariatePoly::new(vec![-x_i, 1.0]))
            .product();

        // denominator
        let denominator = 1.0 / numerator.evaluate(*x);

        numerator.scalar_mul(&denominator)
    }
}

impl Mul for &UnivariatePoly {
    type Output = UnivariatePoly;

    // Multiplying two polynomials or arrays
    fn mul(self, rhs: Self) -> Self::Output {
        // mul for dense
        let new_degree = self.degree() + rhs.degree();
        let mut result = vec![0.0; new_degree + 1];
        for i in 0..self.coefficient.len() {
            for j in 0..rhs.coefficient.len() {
                result[i + j] += self.coefficient[i] * rhs.coefficient[j];
            }
        }
        UnivariatePoly {
            coefficient: result,
        }
    }
}

// This is why it is important to input the array in order of power from smallest to biggest to make it easier to perform actions on them
impl Add for &UnivariatePoly {
    type Output = UnivariatePoly;

    // adding two polynomials
    fn add(self, rhs: Self) -> Self::Output {
        let (mut bigger, smaller) = if self.degree() < rhs.degree() {
            (rhs.clone(), self)
        } else {
            (self.clone(), rhs)
        };

        let _ = bigger
            .coefficient
            .iter_mut()
            .zip(smaller.coefficient.iter())
            .map(|(b_coeff, s_coeff)| *b_coeff += s_coeff)
            .collect::<()>();

        UnivariatePoly::new(bigger.coefficient)
    }
}

impl Sum for UnivariatePoly {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = UnivariatePoly::new(vec![0.0]);
        for item in iter {
            result = &result + &item;
        }
        result
    }
}

impl Product for UnivariatePoly {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = UnivariatePoly::new(vec![1.0]);
        for item in iter {
            result = &result * &item;
        }
        result
    }
}

#[cfg(test)]
mod test {
    use crate::UnivariatePoly;

    fn poly_1() -> UnivariatePoly {
        // f(x) = 1 + 2x + 3x^2
        UnivariatePoly {
            coefficient: vec![1.0, 2.0, 3.0],
        }
    }

    fn poly_2() -> UnivariatePoly {
        // f(x) = 4x + 3 + 5x^11
        UnivariatePoly {
            coefficient: [vec![3.0, 4.0], vec![0.0; 9], vec![5.0]].concat(),
        }
    }

    #[test]
    fn test_degree() {
        assert_eq!(poly_1().degree(), 2);
    }

    #[test]
    fn test_evaluation() {
        assert_eq!(poly_1().evaluate(2.0), 17.0);
    }

    #[test]
    fn test_addition() {
        // f(x) = 1 + 2x + 3x^2
        // f(x) = 4x + 3 + 5x^11

        // r(x) = 4 + 6x + 3x^2 + 5x^11
        assert_eq!(
            (&poly_1() + &poly_2()).coefficient,
            [vec![4.0, 6.0, 3.0], vec![0.0; 8], vec![5.0]].concat()
        )
    }

    #[test]
    fn test_mul() {
        // f(x) = 5 + 2x^2
        let poly_1 = UnivariatePoly {
            coefficient: vec![5.0, 0.0, 2.0],
        };
        // f(x) = 2x + 6
        let poly_2 = UnivariatePoly {
            coefficient: vec![6.0, 2.0],
        };

        // r(x) = 30 + 10x + 12x^2 + 4x^3
        assert_eq!((&poly_1 * &poly_2).coefficient, vec![30.0, 10.0, 12.0, 4.0]);
    }

    #[test]
    fn test_interpolate() {
        // f(x) = 2x
        // [(2, 4), (4, 8)]
        let maybe_2x = UnivariatePoly::interpolate(vec![2.0, 4.0], vec![4.0, 8.0]);
        assert_eq!(maybe_2x.coefficient, vec![0.0, 2.0]);

        // let new_check = UnivariatePoly::interpolate(vec![0.0, 1.0, 2.0, 3.0, 5.0, 10.0], vec![5.0, 7.0, 21.0, 59.0, 255.0, 2005.0]);
        // assert_eq!(new_check.coefficient, vec![5.0, 0.0, 0.0, 2.0]);
    }
}
