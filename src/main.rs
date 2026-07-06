mod bam;
mod barcode;
mod config;
mod count;
mod dedup;
mod filter;
mod gtf;
mod log;
mod orchestrator;
mod trim;

use anyhow::{bail, Result};
use std::time::Instant;

fn main() {
    let run_start = Instant::now();
    match real_main(run_start) {
        Ok(()) => {}
        Err(e) => {
            // Child-process failures already printed a precise FAILURE: line + tail
            // inside checked(); this is the backstop for config/I/O errors. `{e:#}`
            // renders the anyhow chain on a single line — no multi-line Debug dump.
            crate::log::failure("run", &format!("{e:#}"), None);
            std::process::exit(1);
        }
    }
}

fn real_main(start: Instant) -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 || args[1] != "run" || args[2] != "--config" {
        bail!("usage: quartsx run --config <config.yaml>");
    }
    let cfg = config::load(&args[3])?;
    orchestrator::run(&cfg, start)
}
