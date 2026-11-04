#!/usr/bin/env bash

#!/usr/bin/env bash
set -euo pipefail

# Adjust if your cross_sections.xml is in a different subdir:
XS_XML="/xs/endfb-viii.0-hdf5/cross_sections.xml"
CHAIN_XML="/chain/chain_endfb80_fast.xml"

docker run --rm \
  --shm-size=2g \
  -v "$(pwd):/work" \
  -v /opt/data/openmc/xs:/xs:ro \
  -v /opt/data/openmc/chain:/chain:ro \
  -e OPENMC_CROSS_SECTIONS="${XS_XML}" \
  -e OMP_NUM_THREADS=16 \
  -e OPENBLAS_NUM_THREADS=16 \
  -e MKL_NUM_THREADS=16 \
  -w /work \
  prism-openmc \
  bash -lc "test -f '${XS_XML}' && test -f '${CHAIN_XML}' && exec \"\$@\"" -- "$@"
SH

