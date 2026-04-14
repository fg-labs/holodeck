#![deny(unsafe_code)]
#![allow(clippy::cast_precision_loss)]

#[cfg(any(target_pointer_width = "16", target_pointer_width = "32"))]
compile_error!("holodeck requires a 64-bit or wider platform");

pub mod bed;
pub mod commands;
pub mod error_model;
pub mod fasta;
pub mod fragment;
pub mod haplotype;
pub mod output;
pub mod ploidy;
pub mod read;
pub mod read_naming;
pub mod seed;
pub mod sequence_dict;
pub mod vcf;
pub mod version;
