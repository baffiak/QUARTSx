mod bam;
mod barcode;
mod config;
mod count;
mod dedup;
mod filter;
mod gtf;
mod orchestrator;
mod trim;

use anyhow::{bail, Result};

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 || args[1] != "run" || args[2] != "--config" {
        bail!("usage: quartsx run --config <config.yaml>");
    }
    let cfg = config::load(&args[3])?;
    orchestrator::run(&cfg)
}
