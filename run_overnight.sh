#!/usr/bin/env bash
# Overnight batch: matched-k fuel comparison alongside the breakeven equilibrium run.
#
# Usage: nohup ./run_overnight.sh > logs/overnight.log 2>&1 &
#
# Both jobs resume from their last completed depletion step if rerun.  Each job
# is pinned to its own NUMA node (48 threads each, see run_openmc.sh); scaling
# is sublinear, so two jobs give more total throughput than one on all 96.
set -uo pipefail

mkdir -p logs
echo "commit $(git rev-parse HEAD), started $(date -Is)"

(
    for loadout in core_loadouts/matchedk_*.json; do
        NUMA_NODE=1 ./run_openmc.sh \
            python build_prism_core.py "$loadout" --deplete --profile cycle
        echo "$loadout finished (exit $?) $(date -Is)"
    done
    NUMA_NODE=1 ./run_openmc.sh python postprocess.py
) > logs/fuel_comparison.log 2>&1 &
fuel_pid=$!

NUMA_NODE=0 ./run_openmc.sh \
    python -m fuel_cycle.cycle_runner core_loadouts/prism_breakeven_k103.json \
        --profile cycle --max-step-days 170 --k-target 1.03 --k-max 1.035 \
        --cycles 20 --no-stop > logs/breakeven_cycles.log 2>&1
echo "breakeven cycles finished (exit $?) $(date -Is)"

wait "$fuel_pid"
echo "all done $(date -Is)"
