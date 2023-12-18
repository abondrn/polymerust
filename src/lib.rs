#[macro_use]
mod switch_macro;

pub mod polymerust {

use std::cmp;

#[derive(Debug)]
pub struct EnergyParams<T> {
    pub stacking: T,
    pub pair: T,
    pub bulge_loop: T,
    pub internal_loop: T,
}

#[derive(Debug)]
pub struct FoldContext<'a, T> {
    pub sequence: &'a str,
    pub energies: &'a EnergyParams<T>,
    pub mfe: T,
    pub pairs: Vec<(usize, usize)>,
    pub dp: Vec<Vec<T>>,
}

impl<'a, T> FoldContext<'a, T> {
    const EMPTY_VEC: Vec<Vec<T>> = Vec::<Vec::<T>>::new();

    pub fn new(sequence: &'a str, energies: &'a EnergyParams<T>) -> Self {
        Self {
            sequence: sequence,
            energies: energies,
            dp: Self::EMPTY_VEC,
        }
    }
}

impl FoldContext<'_, i32> {

    // Function to calculate RNA secondary structure using Zuker-style dynamic programming
    pub fn zuker(&mut self) {
        let n = self.sequence.len();
        let mut dp = vec![vec![0; n]; n];

        // Fill the DP table using the Zuker-style algorithm
        for len in 1..n {
            for i in 0..(n - len) {
                let j = i + len;
                let mut min_val = i32::MAX;

                let si = self.sequence.as_bytes()[i];
                let sj = self.sequence.as_bytes()[j];

                // Base pairing
                if (si == b'A' && sj == b'U')
                    || (si == b'U' && sj == b'A')
                    || (si == b'C' && sj == b'G')
                    || (si == b'G' && sj == b'C')
                {
                    min_val = cmp::min(min_val, dp[i + 1][j - 1] + self.energies.pair);
                }

                // Hairpin loop
                for k in (i + 3)..j {
                    min_val = cmp::min(min_val, dp[i + 1][k - 1] + dp[k][j] + self.energies.bulge_loop);
                }

                // Internal loop
                if j >= 2 {
                    for k in (i + 1)..(j - 2) {
                        min_val = cmp::min(
                            min_val,
                            dp[i][k] + dp[k + 1][j - 1] + self.energies.internal_loop,
                        );
                    }
                }

                // Stacking
                min_val = cmp::min(min_val, dp[i + 1][j] + self.energies.stacking);
                min_val = cmp::min(min_val, dp[i][j - 1] + self.energies.stacking);

                dp[i][j] = min_val;
            }
        }

        self.dp = dp;
    }

    // Function to traceback and get the secondary structure
    pub fn traceback(&self, i: usize, j: usize) -> String {
        if i >= j {
            return ".".to_string();
        }

        switch! { self.dp[i][j];
            (self.dp[i + 1][j] + self.energies.stacking) => format!(".{}", self.traceback(i + 1, j)),
            (self.dp[i][j - 1] + self.energies.stacking) => format!("{}.", self.traceback(i, j - 1)),
            (self.dp[i + 1][j - 1] + self.energies.pair) => format!("({})", self.traceback(i + 1, j - 1)),
            _ => {

                for k in (i + 3)..j {
                    if self.dp[i][j] == self.dp[i + 1][k - 1] + self.dp[k][j] + self.energies.bulge_loop {
                        return format!(".{}{}", self.traceback(i + 1, k - 1), self.traceback(k, j));
                    }
                }

                for k in (i + 1)..(j - 2) {
                    if self.dp[i][j] == self.dp[i][k] + self.dp[k + 1][j - 1] + self.energies.internal_loop {
                        return format!("{}{}", self.traceback(i, k), self.traceback(k + 1, j - 1));
                    }
                }

                "".to_string()

            },
        }
    }
}

}