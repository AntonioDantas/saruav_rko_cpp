#!/bin/bash
set -e
TIME=${1:-120}
RUNS=${2:-10}
cd "$(dirname "$0")/Program"
make
mkdir -p ../Results
cp config/config_tests.conf config/config_tests.conf.bak
for METHOD in RKO LNS ILS VNS; do
    if [ "$METHOD" = "RKO" ]; then
        LIST="SA ILS VNS BRKGA BRKGA-CS PSO GA LNS GRASP IPR"
    else
        LIST="$METHOD"
    fi
    { printf '%s\n' $LIST; printf '\nMAXRUNS %s\ndebug 0\ncontrol 1\nstrategy 2\nrestart 1\nsizePool 10\n' "$RUNS"; } > config/config_tests.conf
    for f in ../Instances/table2/*.ovrp; do
        ./runTest "table2/$(basename "$f")" "$TIME"
    done
done
mv config/config_tests.conf.bak config/config_tests.conf
