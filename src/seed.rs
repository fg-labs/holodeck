//! Deterministic seed computation for reproducible simulations.
//!
//! When no explicit `--seed` is provided, a seed is derived from the
//! simulation parameters so that identical inputs always produce identical
//! output across runs and processes.

/// Compute a deterministic seed by hashing a string description of the
/// simulation parameters.
///
/// Uses FNV-1a hashing which is deterministic across runs (unlike `ahash`
/// which uses per-process random keys). The same input always produces the
/// same seed.
#[must_use]
pub fn compute_seed(description: &str) -> u64 {
    // FNV-1a hash — deterministic, fast, no external dependencies.
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in description.bytes() {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// Resolve the effective seed: use the explicit seed if provided, otherwise
/// derive one from the given parameter description string.
#[must_use]
pub fn resolve_seed(explicit: Option<u64>, description: &str) -> u64 {
    explicit.unwrap_or_else(|| compute_seed(description))
}

/// Derive a deterministic sub-seed scoped to a named namespace (e.g. a
/// contig name) from a parent seed. Used to produce independent, reproducible
/// RNG streams for per-contig work (reference normalization, etc.) without
/// disturbing the main simulation RNG.
#[must_use]
pub fn derive_seed(parent: u64, namespace: &str) -> u64 {
    compute_seed(&format!("{parent:016x}:{namespace}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_seed_deterministic() {
        let a = compute_seed("hello:42:3.14");
        let b = compute_seed("hello:42:3.14");
        assert_eq!(a, b);
    }

    #[test]
    fn test_compute_seed_known_value() {
        // Pin a known value so we detect if the hash implementation changes.
        let seed = compute_seed("test");
        assert_eq!(seed, 0xf9e6_e6ef_197c_2b25);
    }

    #[test]
    fn test_compute_seed_different_inputs() {
        let a = compute_seed("hello");
        let b = compute_seed("world");
        assert_ne!(a, b);
    }

    #[test]
    fn test_resolve_seed_explicit() {
        let seed = resolve_seed(Some(42), "ignored");
        assert_eq!(seed, 42);
    }

    #[test]
    fn test_resolve_seed_derived() {
        let seed = resolve_seed(None, "hello");
        let expected = compute_seed("hello");
        assert_eq!(seed, expected);
    }

    #[test]
    fn test_derive_seed_deterministic() {
        let a = derive_seed(42, "chr1");
        let b = derive_seed(42, "chr1");
        assert_eq!(a, b);
    }

    #[test]
    fn test_derive_seed_varies_with_namespace() {
        let a = derive_seed(42, "chr1");
        let b = derive_seed(42, "chr2");
        assert_ne!(a, b);
    }

    #[test]
    fn test_derive_seed_varies_with_parent() {
        let a = derive_seed(1, "chr1");
        let b = derive_seed(2, "chr1");
        assert_ne!(a, b);
    }
}
