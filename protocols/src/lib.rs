pub mod fft;
pub mod gkr_2_to_1_trick;
pub mod gkr_circuit;
pub mod gkr_protocol;
pub mod multi_linear;
pub mod partial_sum_check;
pub mod product_poly;
pub mod shamir_secret_sharing;
pub mod sum_check;
pub mod transcript;

pub mod fiat_shamir_non_interactive;
pub mod interactive_sum_check;

use ark_ff::PrimeField;
use std::iter::{Product, Sum};
use std::ops::{Add, Mul};

#[derive(Debug, PartialEq, Clone)]
pub struct UnivariatePoly<F: PrimeField> {
    // 1 coefficient for each power of x
    pub coefficient: Vec<F>,
}

impl<F: PrimeField> UnivariatePoly<F> {
    fn new(coefficient: Vec<F>) -> Self {
        UnivariatePoly { coefficient }
    }

    pub fn degree(&self) -> usize {
        self.coefficient.len() - 1
    }

    pub fn evaluate(&self, x: F) -> F {
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

    pub fn interpolate(xs: &[F], ys: &[F]) -> Self {
        xs.iter()
            .zip(ys.iter())
            .map(|(x, y)| Self::basis(x, &xs).scalar_mul(y))
            .sum()
    }

    // Multiplies an array with an integer
    fn scalar_mul(&self, scalar: &F) -> Self {
        UnivariatePoly {
            coefficient: self.coefficient.iter().map(|x| *x * scalar).collect(),
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
    fn basis(x: &F, interpolating_set: &[F]) -> Self {
        //numerator
        let numerator: UnivariatePoly<F> = interpolating_set
            .iter()
            .filter(|val| *val != x)
            .map(|x_i| UnivariatePoly::new(vec![-*x_i, F::one()]))
            .product();

        // denominator
        let denominator = F::one() / numerator.evaluate(*x);

        numerator.scalar_mul(&denominator)
    }
}

impl<F: PrimeField> Mul for &UnivariatePoly<F> {
    type Output = UnivariatePoly<F>;

    // Multiplying two polynomials or arrays
    fn mul(self, rhs: Self) -> Self::Output {
        // mul for dense
        let new_degree = self.degree() + rhs.degree();
        let mut result = vec![F::zero(); new_degree + 1];
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
impl<F: PrimeField> Add for &UnivariatePoly<F> {
    type Output = UnivariatePoly<F>;

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

impl<F: PrimeField> Sum for UnivariatePoly<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = UnivariatePoly::new(vec![F::zero()]);
        for item in iter {
            result = &result + &item;
        }
        result
    }
}

impl<F: PrimeField> Product for UnivariatePoly<F> {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut result = UnivariatePoly::new(vec![F::one()]);
        for item in iter {
            result = &result * &item;
        }
        result
    }
}

#[cfg(test)]
mod test {
    use crate::UnivariatePoly;
    use ark_bn254::Fq;

    fn poly_1() -> UnivariatePoly<Fq> {
        // f(x) = 1 + 2x + 3x^2
        UnivariatePoly {
            coefficient: vec![Fq::from(1), Fq::from(2), Fq::from(3)],
        }
    }

    fn poly_2() -> UnivariatePoly<Fq> {
        // f(x) = 4x + 3 + 5x^11
        UnivariatePoly {
            coefficient: [
                vec![Fq::from(3), Fq::from(4)],
                vec![Fq::from(0); 9],
                vec![Fq::from(5)],
            ]
            .concat(),
        }
    }

    #[test]
    fn test_degree() {
        assert_eq!(poly_1().degree(), 2);
    }

    #[test]
    fn test_evaluation() {
        assert_eq!(poly_1().evaluate(Fq::from(2)), Fq::from(17));
    }

    #[test]
    fn test_addition() {
        // f(x) = 1 + 2x + 3x^2
        // f(x) = 4x + 3 + 5x^11

        // r(x) = 4 + 6x + 3x^2 + 5x^11
        assert_eq!(
            (&poly_1() + &poly_2()).coefficient,
            [
                vec![Fq::from(4), Fq::from(6), Fq::from(3)],
                vec![Fq::from(0); 8],
                vec![Fq::from(5)]
            ]
            .concat()
        )
    }

    #[test]
    fn test_mul() {
        // f(x) = 5 + 2x^2
        let poly_1 = UnivariatePoly {
            coefficient: vec![Fq::from(5), Fq::from(0), Fq::from(2)],
        };
        // f(x) = 2x + 6
        let poly_2 = UnivariatePoly {
            coefficient: vec![Fq::from(6), Fq::from(2)],
        };

        // r(x) = 30 + 10x + 12x^2 + 4x^3
        assert_eq!(
            (&poly_1 * &poly_2).coefficient,
            vec![Fq::from(30), Fq::from(10), Fq::from(12), Fq::from(4)]
        );
    }

    #[test]
    fn test_interpolate() {
        // f(x) = 2x
        // [(2, 4), (4, 8)]
        let maybe_2x = UnivariatePoly::interpolate(
            &vec![Fq::from(2), Fq::from(4)],
            &vec![Fq::from(4), Fq::from(8)],
        );
        assert_eq!(maybe_2x.coefficient, vec![Fq::from(0), Fq::from(2)]);

        // let new_check = UnivariatePoly::interpolate(vec![0.0, 1.0, 2.0, 3.0, 5.0, 10.0], vec![5.0, 7.0, 21.0, 59.0, 255.0, 2005.0]);
        // assert_eq!(new_check.coefficient, vec![5.0, 0.0, 0.0, 2.0]);
    }

    #[test]
    fn test_fibonacci() {
        // f(x) = 1 + x
        // [(0, 0), (1, 1), (2, 1), (3, 2), (4, 3), (5, 5), (6, 8), (7, 13)]

        let fib = UnivariatePoly::interpolate(
            &vec![
                Fq::from(0),
                Fq::from(1),
                Fq::from(2),
                Fq::from(3),
                Fq::from(4),
                Fq::from(5),
                Fq::from(6),
                Fq::from(7),
            ],
            &vec![
                Fq::from(0),
                Fq::from(1),
                Fq::from(1),
                Fq::from(2),
                Fq::from(3),
                Fq::from(5),
                Fq::from(8),
                Fq::from(13),
            ],
        );

        let check_1 = fib.evaluate(Fq::from(4));
        let check_2 = fib.evaluate(Fq::from(5));
        let check_3 = fib.evaluate(Fq::from(6));
        let check_sum = check_1 + check_2;

        assert_eq!(check_3, check_sum);
    }

    #[test]
    #[should_panic(expected = "assertion `left == right` failed\n  left: 189\n right: 55")]
    fn test_unequal_output() {
        let fib = UnivariatePoly::interpolate(
            &vec![
                Fq::from(0),
                Fq::from(1),
                Fq::from(2),
                Fq::from(3),
                Fq::from(4),
                Fq::from(5),
                Fq::from(6),
                Fq::from(7),
            ],
            &vec![
                Fq::from(0),
                Fq::from(1),
                Fq::from(1),
                Fq::from(2),
                Fq::from(3),
                Fq::from(5),
                Fq::from(8),
                Fq::from(13),
            ],
        );

        let check_1 = fib.evaluate(Fq::from(7));
        let check_2 = fib.evaluate(Fq::from(8));
        let check_3 = fib.evaluate(Fq::from(9));
        let check_sum = check_1 + check_2;

        assert_eq!(check_3, check_sum);
    }

    #[test]
    fn test_gkr_interpolate() {
        let interpolate = UnivariatePoly::interpolate(
            &vec![Fq::from(0), Fq::from(1), Fq::from(2)],
            &vec![Fq::from(0), Fq::from(12), Fq::from(48)],
        );
        dbg!(&interpolate);

        let check_1 = interpolate.evaluate(Fq::from(4));
        dbg!(check_1);
    }
}
