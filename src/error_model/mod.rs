//! Sequencing error models for realistic quality score and error simulation.
//!
//! The [`ErrorModel`] trait defines the interface for position-dependent error
//! models. The [`apply_errors`] free function applies any error model to a
//! read's bases, introducing substitution errors and generating quality scores.

pub mod illumina;

use rand::Rng;

/// Which read in a pair this is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadEnd {
    /// First read in a pair (or the only read for single-end).
    Read1,
    /// Second read in a pair.
    Read2,
}

/// Trait for position-dependent sequencing error models.
///
/// Implementors produce an error probability for each cycle position,
/// allowing realistic simulation of position-dependent error rates.
/// The trait is the extensibility point for future platform support
/// (Element, Ultima, etc.).
pub trait ErrorModel {
    /// Return the probability of a sequencing error at the given cycle
    /// position (0-based) for the given read end.
    fn error_probability(&self, cycle: usize, read_end: ReadEnd) -> f64;

    /// Return the Phred+33 base quality score for the given cycle and read
    /// end (before per-base noise).
    ///
    /// The default implementation computes `Q = -10*log10(p_err)` clamped to
    /// [2, 41] and adds the Phred+33 offset.  Implementors may override this
    /// with a precomputed lookup for speed.
    fn base_quality_phred33(&self, cycle: usize, read_end: ReadEnd) -> u8 {
        let p_err = self.error_probability(cycle, read_end);
        #[expect(clippy::cast_possible_truncation, reason = "clamped to [2, 41]")]
        #[expect(clippy::cast_sign_loss, reason = "clamped to [2, 41]")]
        if p_err > 0.0 {
            (-10.0 * p_err.log10()).round().clamp(2.0, 41.0) as u8 + 33
        } else {
            41 + 33
        }
    }
}

/// Apply sequencing errors to a read's bases in place using the given error
/// model, returning the number of errors and quality scores (Phred+33).
///
/// For each base position:
/// 1. Look up error probability and base quality from the model
/// 2. If a random draw < error probability, substitute a random different base
/// 3. Assign a quality score by adding noise to the base quality
///
/// Returns `(n_errors, quality_scores)`.
pub fn apply_errors(
    model: &impl ErrorModel,
    bases: &mut [u8],
    read_end: ReadEnd,
    rng: &mut impl Rng,
) -> (usize, Vec<u8>) {
    let mut qualities = Vec::with_capacity(bases.len());
    let mut n_errors = 0;

    for (cycle, base) in bases.iter_mut().enumerate() {
        let p_err = model.error_probability(cycle, read_end);

        // Possibly introduce an error.
        if rng.random::<f64>() < p_err {
            *base = random_different_base(*base, rng);
            n_errors += 1;
        }

        // Look up precomputed base quality and add per-base noise.
        let base_qual = model.base_quality_phred33(cycle, read_end);
        let noise: f64 = rng.random_range(-2.0..2.0);
        #[expect(clippy::cast_possible_truncation, reason = "clamped to [35, 74]")]
        #[expect(clippy::cast_sign_loss, reason = "clamped to [35, 74]")]
        let q_noisy = (f64::from(base_qual) + noise).round().clamp(35.0, 74.0) as u8;

        qualities.push(q_noisy);
    }

    (n_errors, qualities)
}

/// Substitute a base with a uniformly random different base.
///
/// Given a DNA base (A, C, G, T), returns one of the other three bases
/// chosen uniformly at random using a single RNG draw and a lookup table.
fn random_different_base(base: u8, rng: &mut impl Rng) -> u8 {
    // For each input base, the three possible substitute bases.
    const ALTS: [u8; 12] = [
        b'C', b'G', b'T', // A -> C, G, T
        b'A', b'G', b'T', // C -> A, G, T
        b'A', b'C', b'T', // G -> A, C, T
        b'A', b'C', b'G', // T -> A, C, G
    ];

    // Map base to row offset: A=0, C=3, G=6, T=9.
    let row = match base {
        b'A' => 0,
        b'C' => 3,
        b'G' => 6,
        _ => 9, // T or any other base
    };

    ALTS[row + rng.random_range(0..3usize)]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_random_different_base_never_same() {
        let mut rng = rand::rng();
        for _ in 0..1000 {
            let original = b'A';
            let result = random_different_base(original, &mut rng);
            assert_ne!(result, original);
            assert!(matches!(result, b'C' | b'G' | b'T'));
        }
    }

    #[test]
    fn test_random_different_base_all_bases() {
        let mut rng = rand::rng();
        for base in [b'A', b'C', b'G', b'T'] {
            let result = random_different_base(base, &mut rng);
            assert_ne!(result, base);
        }
    }
}
