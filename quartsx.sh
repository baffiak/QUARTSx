#!/usr/bin/env bash
# QUARTSx launcher. Usage: quartsx.sh -y config.yaml
#   -y  path to the config yaml
# Run this inside an environment that already provides the dependencies on PATH — the Rust
# toolchain, STAR 2.7.10b, samtools, falco, R and the R packages (e.g. `conda activate quartsx`,
# or any setup you prefer). This script does not create or manage environments.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

yaml=""
while getopts "y:" opt; do
  case "$opt" in
    y) yaml="$OPTARG" ;;
    *) echo "usage: quartsx.sh -y config.yaml" >&2; exit 1 ;;
  esac
done

if [[ -z "$yaml" ]]; then
  echo "usage: quartsx.sh -y config.yaml" >&2
  exit 1
fi
if [[ ! -f "$yaml" ]]; then
  echo "config not found: $yaml" >&2
  exit 1
fi

# build the release binary on first run, using the cargo on PATH
if [[ ! -x "$here/target/release/quartsx" ]]; then
  echo "[quartsx] building release binary" >&2
  cargo build --release --manifest-path "$here/Cargo.toml"
fi

exec "$here/target/release/quartsx" run --config "$yaml"
