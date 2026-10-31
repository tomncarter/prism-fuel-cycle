#!/usr/bin/env bash
# Overnight batch: matched-k fuel comparison alongside the breakeven equilibrium run.
#
# Usage: nohup ./run_overnight.sh > logs/overnight.log 2>&1 &
#
# Both jobs resume from their last completed depletion step if rerun.  Thread
# counts split the 64 hardware threads; scaling is sublinear, so two jobs give
# more total throughput than one job on all of them.
set -uo pipefail

CYCLE_THREADS=44
FUEL_THREADS=20
mkdir -p logs
echo "commit $(git rev-parse HEAD), started $(date -Is)"

(
    for loadout in core_loadouts/matchedk_*.json; do
        ./run_openmc.sh env OMP_NUM_THREADS=$FUEL_THREADS \
            python build_prism_core.py "$loadout" --deplete --profile cycle
        echo "$loadout finished (exit $?) $(date -Is)"
    done
    ./run_openmc.sh python postprocess.py
) > logs/fuel_comparison.log 2>&1 &
fuel_pid=$!

./run_openmc.sh env OMP_NUM_THREADS=$CYCLE_THREADS \
    python -m fuel_cycle.cycle_runner core_loadouts/prism_breakeven_k103.json \
        --profile cycle --max-step-days 170 --k-target 1.03 --k-max 1.035 \
        --cycles 20 --no-stop > logs/breakeven_cycles.log 2>&1
echo "breakeven cycles finished (exit $?) $(date -Is)"

wait "$fuel_pid"
echo "all done $(date -Is)"
