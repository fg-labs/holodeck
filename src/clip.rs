//! Terminal soft-clip artifact model.
//!
//! Real sequencing libraries frequently produce reads whose extreme 5' or 3'
//! ends diverge from the reference enough that an aligner soft-clips them
//! (end-repair fill-in, damaged or otherwise non-templated read ends). Twist
//! EM-seq, for example, shows ~8% of reads carrying a soft clip, strongly
//! 5'-biased and modal at 6-10 bp.
//!
//! Holodeck's simulated reads otherwise align essentially end-to-end, so this
//! module injects that artifact: for a configurable fraction of reads it
//! marks the first/last `K` bases of the *aligned* portion as soft-clipped.
//! The caller corrupts those bases with random sequence (so a downstream
//! aligner re-derives the clip) and records the soft-clip in the golden BAM
//! CIGAR, keeping ground truth exact.
//!
//! The mechanism is protocol-agnostic — it is independent of methylation
//! chemistry and useful for any library type.

use rand::Rng;

/// Configuration for the terminal soft-clip artifact model.
///
/// The 5' and 3' ends clip independently, each with its own probability,
/// because real libraries are markedly asymmetric (the 5' end clips several
/// times more often than the 3' end). The clip-length distribution is shared
/// across both ends, as observed lengths do not differ meaningfully between
/// them.
#[derive(Debug, Clone, Copy)]
pub struct TerminalClipConfig {
    /// Probability in `[0.0, 1.0]` that a read's 5' end receives a clip.
    pub rate_5p: f64,
    /// Probability in `[0.0, 1.0]` that a read's 3' end receives a clip.
    pub rate_3p: f64,
    /// Mean clip length in bases. Lengths are drawn from a truncated
    /// geometric distribution with this mean.
    pub length_mean: usize,
    /// Maximum clip length in bases; sampled lengths are clamped to this.
    pub length_max: usize,
}

impl TerminalClipConfig {
    /// Whether either end can produce a clip. When this is `false`,
    /// [`Self::sample_clips`] consumes no randomness and returns `(0, 0)`, so
    /// the simulator's output is byte-identical to having no clip model at
    /// all.
    #[must_use]
    pub fn is_enabled(&self) -> bool {
        self.rate_5p > 0.0 || self.rate_3p > 0.0
    }

    /// Sample `(clip_5p, clip_3p)` lengths for one read whose aligned
    /// (genomic) portion is `aligned_len` bases long.
    ///
    /// Each end independently clips with its configured probability; a chosen
    /// clip's length is drawn from a truncated geometric with mean
    /// `length_mean`, clamped to `[1, length_max]`. The combined clip is then
    /// capped so at least one aligned base survives (an alignment needs a
    /// match block), favouring the 5' end when both are drawn and the budget
    /// is tight.
    ///
    /// Returns `(0, 0)` without drawing any randomness when the model is
    /// disabled or `aligned_len < 2`.
    pub fn sample_clips(&self, aligned_len: usize, rng: &mut impl Rng) -> (usize, usize) {
        // No clip is representable without leaving an aligned base behind, and
        // a disabled model must not perturb the RNG stream.
        if !self.is_enabled() || aligned_len < 2 {
            return (0, 0);
        }

        let mut clip_5p = self.sample_end(self.rate_5p, rng);
        let mut clip_3p = self.sample_end(self.rate_3p, rng);

        // Cap so at least one aligned base remains. Trim the 3' end first
        // (the rarer, less characteristic artifact), then the 5' end.
        let budget = aligned_len - 1;
        if clip_5p + clip_3p > budget {
            clip_3p = clip_3p.min(budget.saturating_sub(clip_5p));
            clip_5p = clip_5p.min(budget);
        }

        (clip_5p, clip_3p)
    }

    /// Draw a clip length for a single end: `0` if the end does not clip,
    /// otherwise a truncated-geometric length in `[1, length_max]`. Consumes
    /// no randomness when `rate <= 0.0`.
    fn sample_end(&self, rate: f64, rng: &mut impl Rng) -> usize {
        if rate <= 0.0 {
            return 0;
        }
        if rng.random::<f64>() >= rate {
            return 0;
        }
        self.sample_length(rng)
    }

    /// Sample a clip length from a truncated geometric distribution with mean
    /// [`Self::length_mean`], clamped to `[1, length_max]`.
    fn sample_length(&self, rng: &mut impl Rng) -> usize {
        let max = self.length_max.max(1);
        let mean = self.length_mean.max(1);
        if mean <= 1 {
            // Degenerate distribution: always the shortest representable clip.
            return 1;
        }
        #[expect(clippy::cast_precision_loss, reason = "small clip means")]
        let p = 1.0 / mean as f64;
        let u: f64 = rng.random(); // [0, 1)
        // Inverse-CDF of the geometric (trials >= 1): k = floor(ln(1-u)/ln(1-p)) + 1.
        let k = ((1.0 - u).ln() / (1.0 - p).ln()).floor();
        #[expect(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "k is finite, non-negative after floor"
        )]
        let k = k as usize + 1;
        k.clamp(1, max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    fn cfg(rate_5p: f64, rate_3p: f64) -> TerminalClipConfig {
        TerminalClipConfig { rate_5p, rate_3p, length_mean: 8, length_max: 20 }
    }

    #[test]
    fn disabled_config_is_not_enabled() {
        assert!(!cfg(0.0, 0.0).is_enabled());
        assert!(cfg(0.1, 0.0).is_enabled());
        assert!(cfg(0.0, 0.1).is_enabled());
    }

    #[test]
    fn disabled_config_returns_zero_and_consumes_no_randomness() {
        // A disabled model must leave the RNG stream untouched so existing
        // simulations stay byte-identical.
        let mut rng_a = SmallRng::seed_from_u64(7);
        let clips = cfg(0.0, 0.0).sample_clips(100, &mut rng_a);
        assert_eq!(clips, (0, 0));
        let after: u64 = rng_a.random();

        let mut rng_b = SmallRng::seed_from_u64(7);
        let direct: u64 = rng_b.random();
        assert_eq!(after, direct, "disabled clip model must not advance the RNG");
    }

    #[test]
    fn short_alignment_cannot_clip() {
        let mut rng = SmallRng::seed_from_u64(1);
        assert_eq!(cfg(1.0, 1.0).sample_clips(0, &mut rng), (0, 0));
        assert_eq!(cfg(1.0, 1.0).sample_clips(1, &mut rng), (0, 0));
    }

    #[test]
    fn certain_5p_clip_always_clips_5p_only() {
        let mut rng = SmallRng::seed_from_u64(42);
        for _ in 0..200 {
            let (c5, c3) = cfg(1.0, 0.0).sample_clips(150, &mut rng);
            assert!(c5 >= 1, "rate_5p=1.0 must always clip the 5' end");
            assert_eq!(c3, 0, "rate_3p=0.0 must never clip the 3' end");
            assert!(c5 <= 20, "clip must respect length_max");
        }
    }

    #[test]
    fn certain_3p_clip_always_clips_3p_only() {
        let mut rng = SmallRng::seed_from_u64(42);
        for _ in 0..200 {
            let (c5, c3) = cfg(0.0, 1.0).sample_clips(150, &mut rng);
            assert_eq!(c5, 0);
            assert!((1..=20).contains(&c3));
        }
    }

    #[test]
    fn clip_never_consumes_whole_alignment() {
        let mut rng = SmallRng::seed_from_u64(99);
        // length_max larger than the alignment forces the cap to engage.
        let config =
            TerminalClipConfig { rate_5p: 1.0, rate_3p: 1.0, length_mean: 8, length_max: 50 };
        for aligned in 2..30 {
            for _ in 0..50 {
                let (c5, c3) = config.sample_clips(aligned, &mut rng);
                assert!(c5 + c3 <= aligned - 1, "must leave >=1 aligned base (aligned={aligned})");
            }
        }
    }

    #[test]
    fn observed_rate_tracks_configured_rate() {
        let mut rng = SmallRng::seed_from_u64(2024);
        let mut clipped = 0;
        let n = 20_000;
        for _ in 0..n {
            let (c5, _) = cfg(0.3, 0.0).sample_clips(150, &mut rng);
            if c5 > 0 {
                clipped += 1;
            }
        }
        let observed = f64::from(clipped) / f64::from(n);
        assert!((observed - 0.3).abs() < 0.02, "expected ~0.30 clip rate, got {observed}");
    }

    #[test]
    fn mean_length_is_in_the_right_ballpark() {
        let mut rng = SmallRng::seed_from_u64(2025);
        let config =
            TerminalClipConfig { rate_5p: 1.0, rate_3p: 0.0, length_mean: 8, length_max: 1000 };
        let mut sum = 0usize;
        let n = 50_000;
        for _ in 0..n {
            let (c5, _) = config.sample_clips(2000, &mut rng);
            sum += c5;
        }
        #[expect(clippy::cast_precision_loss, reason = "test arithmetic")]
        let mean = sum as f64 / f64::from(n);
        assert!((mean - 8.0).abs() < 1.0, "expected mean ~8, got {mean}");
    }
}
