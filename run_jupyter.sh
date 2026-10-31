#!/usr/bin/env bash
set -euo pipefail

PORT="${1:-8888}"

docker run --rm \
  -p "${PORT}:8888" \
  -v "$(pwd):/work" \
  -w /work \
  prism-openmc \
  jupyter lab \
    --ip=0.0.0.0 \
    --port=8888 \
    --no-browser \
    --notebook-dir=/work \
    --ServerApp.token='' \
    --ServerApp.password=''
