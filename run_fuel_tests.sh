#!/usr/bin/env bash
# Run every core_loadouts/fueltest_*.json at a given profile, then postprocess.
#
# Usage:
#   ./run_fuel_tests.sh [profile] [--deplete]
#
#   ./run_fuel_tests.sh smoke              # pipeline check, transport only
#   ./run_fuel_tests.sh quick              # BOC k-eff screen for all fuels
#   ./run_fuel_tests.sh quick --deplete    # ~1-year depletion, 6 steps each
set -uo pipefail

PROFILE="${1:-quick}"
shift || true

failed=()
for loadout in core_loadouts/fueltest_*.json; do
    ./run_openmc.sh python build_prism_core.py "$loadout" --profile "$PROFILE" "$@"
    rc=$?
    # 127 comes from the stray trailing line in run_openmc.sh and is harmless.
    if [[ $rc -ne 0 && $rc -ne 127 ]]; then
        echo "!! $loadout failed (exit $rc)"
        failed+=("$loadout")
    fi
done

./run_openmc.sh python postprocess.py

if (( ${#failed[@]} )); then
    echo "Failed loadouts: ${failed[*]}"
    exit 1
fi
