use super::{ErrorModel, ReadEnd};

/// Fraction of read length at which the quality decay ramp begins.
const DECAY_START_FRACTION: f64 = 0.7;

/// Default R2 error rate multiplier relative to R1.
const DEFAULT_R2_MULTIPLIER: f64 = 1.5;

/// Precomputed per-cycle error probability and base quality score.
#[derive(Debug, Clone, Copy)]
struct CycleParams {
    /// Error probability at this cycle.
    p_err: f64,
    /// Base quality score (before noise): Q = -10*log10(p_err), Phred+33.
    base_qual: u8,
}

/// Illumina-style position-dependent error model.
///
/// Error rate starts at `min_error_rate` for early cycles and linearly ramps
/// to `max_error_rate` by the last cycle. The ramp begins at approximately
/// 70% of the read length, mimicking the quality decay caused by phasing and
/// dephasing in Illumina sequencing.
///
/// Read 2 has a higher error rate than read 1 (controlled by
/// `r2_multiplier`), reflecting the typically worse quality of the second
/// read in a pair.
///
/// Error probabilities and base quality scores are precomputed at
/// construction time for each cycle and read end, so the per-base hot path
/// is a table lookup rather than repeated floating-point math.
pub struct IlluminaErrorModel {
    /// Read length this model is configured for.
    read_length: usize,
    /// Error rate at the start of reads (e.g. 0.001).
    min_error_rate: f64,
    /// Maximum error rate at the end of reads (e.g. 0.01).
    max_error_rate: f64,
    /// Multiplier applied to the entire R2 error curve (e.g. 1.5).
    r2_multiplier: f64,
    /// Cycle position where the error ramp begins (computed from read_length).
    decay_start: usize,
    /// Precomputed per-cycle parameters for R1.
    r1_params: Vec<CycleParams>,
    /// Precomputed per-cycle parameters for R2.
    r2_params: Vec<CycleParams>,
}

impl IlluminaErrorModel {
    /// Create a new Illumina error model with default R2 multiplier (1.5x).
    ///
    /// # Arguments
    /// * `read_length` — Number of cycles (bases) per read.
    /// * `min_error_rate` — Per-base error rate at the start of reads.
    /// * `max_error_rate` — Per-base error rate at the end of reads.
    ///
    /// # Panics
    /// Panics if `read_length` is 0 or if `min_error_rate > max_error_rate`.
    #[must_use]
    pub fn new(read_length: usize, min_error_rate: f64, max_error_rate: f64) -> Self {
        Self::with_r2_multiplier(read_length, min_error_rate, max_error_rate, DEFAULT_R2_MULTIPLIER)
    }

    /// Create a new Illumina error model with a custom R2 multiplier.
    ///
    /// The `r2_multiplier` scales the entire error curve for read 2 relative
    /// to read 1. A value of 1.0 means R1 and R2 have identical error rates.
    ///
    /// # Panics
    /// Panics if `read_length` is 0 or if `min_error_rate > max_error_rate`.
    #[must_use]
    pub fn with_r2_multiplier(
        read_length: usize,
        min_error_rate: f64,
        max_error_rate: f64,
        r2_multiplier: f64,
    ) -> Self {
        assert!(read_length > 0, "read_length must be > 0");
        assert!(
            min_error_rate <= max_error_rate,
            "min_error_rate ({min_error_rate}) must be <= max_error_rate ({max_error_rate})"
        );

        #[expect(clippy::cast_possible_truncation, reason = "product is bounded by read_length")]
        #[expect(clippy::cast_sign_loss, reason = "product is non-negative")]
        let decay_start = (read_length as f64 * DECAY_START_FRACTION).round() as usize;

        let mut model = Self {
            read_length,
            min_error_rate,
            max_error_rate,
            r2_multiplier,
            decay_start,
            r1_params: Vec::new(),
            r2_params: Vec::new(),
        };

        // Precompute per-cycle error probabilities and base quality scores.
        model.r1_params =
            (0..read_length).map(|c| model.compute_cycle_params(c, ReadEnd::Read1)).collect();
        model.r2_params =
            (0..read_length).map(|c| model.compute_cycle_params(c, ReadEnd::Read2)).collect();

        model
    }
}

impl IlluminaErrorModel {
    /// Compute error probability for a given cycle using the ramp model.
    fn compute_error_probability(&self, cycle: usize, read_end: ReadEnd) -> f64 {
        let base_rate = if cycle < self.decay_start {
            self.min_error_rate
        } else {
            let ramp_len = self.read_length.saturating_sub(self.decay_start).max(1);
            let progress = (cycle - self.decay_start) as f64 / ramp_len as f64;
            self.min_error_rate + (self.max_error_rate - self.min_error_rate) * progress.min(1.0)
        };

        match read_end {
            ReadEnd::Read1 => base_rate,
            ReadEnd::Read2 => (base_rate * self.r2_multiplier).min(1.0),
        }
    }

    /// Precompute the error probability and base quality for a cycle.
    fn compute_cycle_params(&self, cycle: usize, read_end: ReadEnd) -> CycleParams {
        let p_err = self.compute_error_probability(cycle, read_end);

        #[expect(clippy::cast_possible_truncation, reason = "clamped to [2, 41]")]
        #[expect(clippy::cast_sign_loss, reason = "clamped to [2, 41]")]
        let base_qual = if p_err > 0.0 {
            (-10.0 * p_err.log10()).round().clamp(2.0, 41.0) as u8 + 33
        } else {
            41 + 33
        };

        CycleParams { p_err, base_qual }
    }

    /// Return the precomputed cycle parameters for the given cycle and read end.
    #[inline]
    fn cycle_params(&self, cycle: usize, read_end: ReadEnd) -> CycleParams {
        match read_end {
            ReadEnd::Read1 => self.r1_params[cycle],
            ReadEnd::Read2 => self.r2_params[cycle],
        }
    }
}

impl ErrorModel for IlluminaErrorModel {
    fn error_probability(&self, cycle: usize, read_end: ReadEnd) -> f64 {
        self.cycle_params(cycle, read_end).p_err
    }

    fn base_quality_phred33(&self, cycle: usize, read_end: ReadEnd) -> u8 {
        self.cycle_params(cycle, read_end).base_qual
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_rate_at_start() {
        let model = IlluminaErrorModel::new(150, 0.001, 0.01);
        let p = model.error_probability(0, ReadEnd::Read1);
        assert!((p - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_error_rate_flat_before_decay() {
        let model = IlluminaErrorModel::new(150, 0.001, 0.01);
        // decay_start ≈ 105 for read_length 150
        for cycle in 0..100 {
            let p = model.error_probability(cycle, ReadEnd::Read1);
            assert!((p - 0.001).abs() < 1e-10, "cycle {cycle}: expected 0.001, got {p}");
        }
    }

    #[test]
    fn test_error_rate_ramps_to_max() {
        let model = IlluminaErrorModel::new(150, 0.001, 0.01);
        let p = model.error_probability(149, ReadEnd::Read1);
        assert!((p - 0.01).abs() < 0.001, "expected ~0.01 at last cycle, got {p}");
    }

    #[test]
    fn test_error_rate_increases_monotonically() {
        let model = IlluminaErrorModel::new(150, 0.001, 0.01);
        let mut prev = 0.0;
        for cycle in 0..150 {
            let p = model.error_probability(cycle, ReadEnd::Read1);
            assert!(p >= prev, "cycle {cycle}: {p} < {prev}");
            prev = p;
        }
    }

    #[test]
    fn test_r2_higher_than_r1() {
        let model = IlluminaErrorModel::new(150, 0.001, 0.01);
        for cycle in [0, 50, 100, 149] {
            let p1 = model.error_probability(cycle, ReadEnd::Read1);
            let p2 = model.error_probability(cycle, ReadEnd::Read2);
            assert!(p2 > p1, "cycle {cycle}: R2 ({p2}) should be > R1 ({p1})");
        }
    }

    #[test]
    fn test_r2_multiplier() {
        let model = IlluminaErrorModel::new(150, 0.001, 0.01);
        let p1 = model.error_probability(0, ReadEnd::Read1);
        let p2 = model.error_probability(0, ReadEnd::Read2);
        assert!((p2 - p1 * 1.5).abs() < 1e-10, "R2 should be 1.5x R1 at start");
    }

    #[test]
    fn test_apply_errors_reproducible() {
        use rand::SeedableRng;
        use rand::rngs::SmallRng;

        let model = IlluminaErrorModel::new(10, 0.1, 0.3);
        let mut bases1 = b"ACGTACGTAC".to_vec();
        let mut bases2 = b"ACGTACGTAC".to_vec();

        let mut rng1 = SmallRng::seed_from_u64(42);
        let mut rng2 = SmallRng::seed_from_u64(42);

        let (n1, q1) =
            crate::error_model::apply_errors(&model, &mut bases1, ReadEnd::Read1, &mut rng1);
        let (n2, q2) =
            crate::error_model::apply_errors(&model, &mut bases2, ReadEnd::Read1, &mut rng2);

        assert_eq!(bases1, bases2);
        assert_eq!(q1, q2);
        assert_eq!(n1, n2);
    }

    #[test]
    fn test_apply_errors_introduces_errors() {
        use rand::SeedableRng;
        use rand::rngs::SmallRng;

        let model = IlluminaErrorModel::new(100, 0.5, 0.5); // High error rate
        let original = vec![b'A'; 100];
        let mut bases = original.clone();
        let mut rng = SmallRng::seed_from_u64(123);

        let (n_errors, qualities) =
            crate::error_model::apply_errors(&model, &mut bases, ReadEnd::Read1, &mut rng);

        assert!(n_errors > 10, "expected many errors with 50% rate, got {n_errors}");
        assert_eq!(qualities.len(), 100);
        // Some bases should have changed.
        assert_ne!(bases, original);
    }

    #[test]
    fn test_apply_errors_quality_encoding() {
        use rand::SeedableRng;
        use rand::rngs::SmallRng;

        let model = IlluminaErrorModel::new(10, 0.001, 0.01);
        let mut bases = b"ACGTACGTAC".to_vec();
        let mut rng = SmallRng::seed_from_u64(99);

        let (_, qualities) =
            crate::error_model::apply_errors(&model, &mut bases, ReadEnd::Read1, &mut rng);

        for &q in &qualities {
            // Phred+33: quality 2..=41 maps to ASCII 35..=74
            assert!(q >= 35, "quality {q} below minimum Phred+33 (35)");
            assert!(q <= 74, "quality {q} above maximum Phred+33 (74)");
        }
    }

    #[test]
    fn test_apply_errors_zero_rate() {
        use rand::SeedableRng;
        use rand::rngs::SmallRng;

        let model = IlluminaErrorModel::new(10, 0.0, 0.0);
        let original = b"ACGTACGTAC".to_vec();
        let mut bases = original.clone();
        let mut rng = SmallRng::seed_from_u64(42);

        let (n_errors, qualities) =
            crate::error_model::apply_errors(&model, &mut bases, ReadEnd::Read1, &mut rng);

        assert_eq!(n_errors, 0, "zero error rate should produce no errors");
        assert_eq!(bases, original, "bases should be unchanged");
        assert_eq!(qualities.len(), 10);
    }

    #[test]
    #[should_panic(expected = "read_length must be > 0")]
    fn test_zero_read_length_panics() {
        let _ = IlluminaErrorModel::new(0, 0.001, 0.01);
    }

    #[test]
    #[should_panic(expected = "min_error_rate")]
    fn test_min_gt_max_panics() {
        let _ = IlluminaErrorModel::new(150, 0.1, 0.01);
    }
}
