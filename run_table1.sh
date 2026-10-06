#!/bin/bash
set -e
TIME=${1:-120}
RUNS=${2:-10}
cd "$(dirname "$0")/Program"
make
mkdir -p ../Results
cp config/config_tests.conf config/config_tests.conf.bak
{ printf '%s\n' SA ILS VNS BRKGA BRKGA-CS PSO GA LNS GRASP IPR; printf '\nMAXRUNS %s\ndebug 0\ncontrol 1\nstrategy 2\nrestart 1\nsizePool 10\n' "$RUNS"; } > config/config_tests.conf
for f in ../Instances/table1/*.ovrp; do
    ./runTest "table1/$(basename "$f")" "$TIME"
done
mv config/config_tests.conf.bak config/config_tests.conf
