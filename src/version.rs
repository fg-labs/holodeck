#![allow(clippy::doc_markdown)] // Generated file contains OPT_LEVEL without backticks

use std::sync::LazyLock;

include!(concat!(env!("OUT_DIR"), "/built.rs"));

/// Version of the software including:
/// - Cargo package version
/// - Short git commit hash (if available)
/// - Git dirty flag (if the repo had uncommitted changes at build time)
pub static VERSION: LazyLock<String> = LazyLock::new(|| {
    let prefix = if let Some(hash) = GIT_COMMIT_HASH {
        let short = &hash[..hash.len().min(8)];
        format!("{PKG_VERSION}-{short}")
    } else {
        PKG_VERSION.to_string()
    };
    let suffix = match GIT_DIRTY {
        Some(true) => "-dirty",
        _ => "",
    };
    format!("{prefix}{suffix}")
});
