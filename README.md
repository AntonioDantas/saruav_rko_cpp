# RKO for Multi-UAV Search and Rescue

Code and instances for **"Maximizing Detection Probability in UAV Search Missions using Random-Key Optimizers"**.

The solver is the Random-Key Optimizer (RKO) framework v2.0 in C++, with the reader and decoder of the probabilistic multi-depot Open Vehicle Routing Problem studied in the article.

## Contents

```
Program/
  src/Main/main.cpp            RKO driver
  src/Method.h                 shared search components
  src/QLearning.h              online parameter control
  src/MH/                      BRKGA, BRKGA-CS, SA, ILS, VNS, LNS, PSO, GA, GRASP, IPR, MultiStart
  src/Problem/Problem.h        instance reader and decoder
  src/Output.h                 result files
  config/config_tests.conf     metaheuristics, runs and search settings
  config/Parameters*.txt       parameter sets
Instances/
  table1/                      exact method x RKO comparison (10, 50 and 100 POIs)
  table2/                      RKO x single metaheuristics comparison (100 POIs)
run_table1.sh
run_table2.sh
```

## Requirements

A C++20 compiler with OpenMP.

- Linux: `g++`
- macOS: GCC from Homebrew (`brew install gcc`); Apple clang has no OpenMP
- Windows: MinGW-w64 `g++`

## Build

```bash
cd Program
make
```

On macOS, `make` uses `g++-14`; for another version pass it explicitly, e.g. `make CCC=g++-15`.

Without `make`:

```bash
g++ -fopenmp -O3 -std=c++2a -o runTest src/Main/main.cpp
```

## Run

From `Program/`, give the instance path relative to `Instances/` and the time limit in seconds:

```bash
./runTest "table1/1000N_10p_0to3_3d_60E_f1=40_f2=0.ovrp" 120
```

To reproduce the tables (120 s per run, 10 runs per instance):

```bash
./run_table1.sh
./run_table2.sh
```

Both accept the time limit and the number of runs: `./run_table1.sh 60 3`.

`run_table1.sh` runs the full RKO. `run_table2.sh` runs the full RKO and then LNS, ILS and VNS alone, each on one thread.

## Output

Results are appended to `Results/`:

- `Results_RKO.csv`, one line per instance and method: metaheuristics used, best objective value, its f1 (normalized travel) and f2 (collected probability), mean objective, time to best, total time and the objective of every run.
- `Solutions_RKO.txt`: best random-key vector and a per-UAV summary.

The objective is minimized. For Model 1 instances it is the negative of the normalized collected probability `Z1`; for Model 2 it is the normalized travel `Z2`.

## Configuration

`config/config_tests.conf`:

| Entry | Meaning |
|---|---|
| metaheuristic names | one per line, each runs on its own thread and shares the elite pool |
| `MAXRUNS` | independent runs per instance |
| `debug` | 0 writes the result files with a clock seed; 1 prints to the screen with a fixed seed |
| `control` | 0 offline parameters, 1 online Q-learning control |
| `strategy` | local search: 1 first improvement, 2 best improvement |
| `restart` | fraction of the time limit after which the search restarts |
| `sizePool` | size of the elite pool |

## Instance format

```
Nodes 13                 depots + POIs
Vehicles 3               UAVs, one depot each
Vmin 0                   minimum visits per POI
Vmax 3                   maximum visits per POI
F1max 40                 T_max per UAV (Model 1); 0 selects Model 2
F2max 0                  P_min (Model 2)
Id,Probabilties,Energy
0 0.0 37                 depot: energy > 0 is the battery E_k of its UAV
3 0.0604 0.0             POI: probability p_i
...
Distances                travel-cost matrix between all nodes
```

The sections after the matrix (grid cells, coordinates and corridor paths) describe the geography and are not read by the solver. File names follow `1000N_<POIs>p_<Vmin>to<Vmax>_<UAVs>d_<E>E_f1=<T_max>_f2=<P_min>`.

The `table2` instances correspond to the article as follows:

| Article | File |
|---|---|
| i04 | `1000N_100p_0to1_5d_500E_f1=100_f2=0` |
| i05 | `1000N_100p_1to2_5d_500E_f1=300_f2=0` |
| i06 | `1000N_100p_0to1_3d_500E_f1=100_f2=0` |
| i07 | `1000N_100p_1to2_3d_500E_f1=300_f2=0` |
| i08 | `1000N_100p_0to1_5d_500E_f1=0_f2=0.8` |
| i09 | `1000N_100p_1to2_5d_500E_f1=0_f2=0.8` |
| i10 | `1000N_100p_0to1_3d_500E_f1=0_f2=0.8` |
| i11 | `1000N_100p_1to2_3d_500E_f1=0_f2=0.8` |

## License

MIT, see `LICENSE`. When using the RKO, please cite:

Chaves, A.A., Resende, M.G.C., Schuetz, M.J.A., Brubaker, J.K., Katzgraber, H.G., Arruda, E.F., Silva, R.M.A. A Random-Key Optimizer for Combinatorial Optimization. *Journal of Heuristics* 31(4), 32, 2025. https://doi.org/10.1007/s10732-025-09568-z
