#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$project_dir/logs"
Rscript "$project_dir/workflow/10_cellchat.R" "$@" 2>&1 | tee "$project_dir/logs/10_cellchat.log"
