#![deny(unsafe_code)]
#![allow(clippy::cast_precision_loss)]

use std::io::Write;
use std::time::Instant;

use anyhow::Result;
use chrono::Local;
use clap::Parser;
use clap::builder::styling::{AnsiColor, Effects, Styles};
use env_logger::Env;
use holodeck_lib::commands::command::Command;
use holodeck_lib::commands::eval::Eval;
use holodeck_lib::commands::methylate::Methylate;
use holodeck_lib::commands::mutate::Mutate;
use holodeck_lib::commands::simulate::Simulate;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Cyan.on_default())
    .valid(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .invalid(AnsiColor::Yellow.on_default().effects(Effects::BOLD));

/// Modern NGS read simulator.
///
/// Holodeck generates simulated sequencing reads from a reference genome with
/// optional variants from a VCF file.  It supports Illumina-style paired-end and
/// single-end reads with position-dependent error profiles.
#[derive(Parser)]
#[command(
    name = "holodeck",
    version = holodeck_lib::version::VERSION.as_str(),
    long_about,
    styles = STYLES,
    after_long_help = "Run 'holodeck <command> --help' for detailed options on each command."
)]
struct Cli {
    /// Enable verbose (debug-level) logging.
    #[arg(long, global = true)]
    verbose: bool,

    #[command(subcommand)]
    command: Subcommand,
}

/// All holodeck subcommands.
#[derive(clap::Subcommand)]
enum Subcommand {
    Simulate(Simulate),
    Mutate(Mutate),
    Methylate(Methylate),
    Eval(Eval),
}

impl Command for Subcommand {
    fn execute(&self) -> Result<()> {
        match self {
            Subcommand::Simulate(c) => c.execute(),
            Subcommand::Mutate(c) => c.execute(),
            Subcommand::Methylate(c) => c.execute(),
            Subcommand::Eval(c) => c.execute(),
        }
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let log_level = if cli.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(Env::default().default_filter_or(log_level))
        .format(|buf, record| {
            let timestamp = Local::now().format("%Y-%m-%d %H:%M:%S");
            writeln!(buf, "[{timestamp} {} {}] {}", record.level(), record.target(), record.args())
        })
        .init();

    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    log::info!("Holodeck by Fulcrum Genomics - https://www.github.com/fulcrumgenomics/holodeck");
    log::info!("Executing: {cmdline}");

    let start = Instant::now();
    let result = cli.command.execute();
    let elapsed = start.elapsed();
    let minutes = elapsed.as_secs() / 60;
    let seconds = elapsed.as_secs() % 60;

    match &result {
        Ok(()) => log::info!("Successfully completed execution in {minutes}m:{seconds:02}s."),
        Err(_) => log::error!("Execution failed after {minutes}m:{seconds:02}s."),
    }

    result
}
