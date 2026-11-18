#!/usr/bin/env bash
set -euo pipefail

# Adjust if your cross_sections.xml is in a different subdir:
XS_XML="/xs/endfb-viii.0-hdf5/cross_sections.xml"
CHAIN_XML="/chain/chain_endfb80_fast.xml"

# Pin each run to one NUMA node (CPUs + memory) so two runs can go side by side:
#   NUMA_NODE=0 ./run_openmc.sh ...   # default
#   NUMA_NODE=1 ./run_openmc.sh ...
# THREADS defaults to every logical CPU on that node.
NUMA_NODE="${NUMA_NODE:-0}"
NODE_CPUS="$(cat "/sys/devices/system/node/node${NUMA_NODE}/cpulist")"
NODE_NCPU="$(lscpu -p=CPU,NODE | grep -v '^#' | awk -F, -v n="${NUMA_NODE}" '$2==n' | wc -l)"
THREADS="${THREADS:-${NODE_NCPU}}"

docker run --rm \
  --shm-size=2g \
  --cpuset-cpus="${NODE_CPUS}" \
  --cpuset-mems="${NUMA_NODE}" \
  -v "$(pwd):/work" \
  -v /opt/data/openmc/xs:/xs:ro \
  -v /opt/data/openmc/chain:/chain:ro \
  -e OPENMC_CROSS_SECTIONS="${XS_XML}" \
  -e OMP_NUM_THREADS="${THREADS}" \
  -e OPENBLAS_NUM_THREADS="${THREADS}" \
  -e MKL_NUM_THREADS="${THREADS}" \
  -w /work \
  prism-openmc \
  bash -lc "test -f '${XS_XML}' && test -f '${CHAIN_XML}' && exec \"\$@\"" -- "$@"
