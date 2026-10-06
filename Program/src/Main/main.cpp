#include <sys/time.h>
#include <math.h>
#include <cstring>
#include <ctime>
#include <iostream>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>

#else
#include <time.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <algorithm>
#include <utility>
#include <numeric>
#include <map>
#include <limits>
#include <random>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <atomic>

#include "../Data.h"

std::mt19937 rng;

std::atomic<bool> stop_execution(false);

std::vector<TSol> pool;

#include "../Problem/Problem.h"
#include "../Output.h"
#include "../Method.h"
#include "../QLearning.h"
#include "../MH/BRKGA.h"
#include "../MH/SA.h"
#include "../MH/GRASP.h"
#include "../MH/ILS.h"
#include "../MH/VNS.h"
#include "../MH/PSO.h"
#include "../MH/GA.h"
#include "../MH/BRKGA_CS.h"
#include "../MH/LNS.h"
#include "../MH/IPR.h"
#include "../MH/MultiStart.h"

int main(int argc, char *argv[])
{

    TRunData runData;

    char nameInstance[256];

    strncpy(nameInstance, argv[1], 255);
    runData.MAXTIME = std::stoi(argv[2]);

#define TOTAL_MH 11

    const char *all_algorithms[TOTAL_MH] = {
        "BRKGA", "SA", "GRASP", "ILS", "VNS", "PSO", "GA", "LNS", "BRKGA-CS", "MultiStart", "IPR"};

    void BRKGA(const TRunData &, const TProblemData &);
    void SA(const TRunData &, const TProblemData &);
    void GRASP(const TRunData &, const TProblemData &);
    void ILS(const TRunData &, const TProblemData &);
    void VNS(const TRunData &, const TProblemData &);
    void PSO(const TRunData &, const TProblemData &);
    void GA(const TRunData &, const TProblemData &);
    void LNS(const TRunData &, const TProblemData &);
    void BRKGA_CS(const TRunData &, const TProblemData &);
    void MultiStart(const TRunData &, const TProblemData &);
    void IPR(const TRunData &, const TProblemData &);

    void (*all_functions[TOTAL_MH])(const TRunData &, const TProblemData &) = {
        BRKGA, SA, GRASP, ILS, VNS, PSO, GA, LNS, BRKGA_CS, MultiStart, IPR};

#define MAX_MH 100
    void (*functions_MH[MAX_MH])(const TRunData &, const TProblemData &);
    const char *algorithms[MAX_MH];
    int NUM_MH = 0;

    FILE *fileConf = fopen("config/config_tests.conf", "r");
    if (!fileConf)
    {
        printf("\nERROR: File config-tests.conf not found\n");
        exit(1);
    }

    char line[250];
    runData.MAXRUNS = 0;
    runData.debug = 0;
    runData.control = 0;

    while (fgets(line, sizeof(line), fileConf))
    {

        line[strcspn(line, "\n")] = '\0';

        int idx = -1;
        for (int i = 0; i < TOTAL_MH; ++i)
        {
            if (strcmp(line, all_algorithms[i]) == 0)
            {
                idx = i;
            }
        }

        if (idx != -1)
        {
            functions_MH[NUM_MH] = all_functions[idx];
            algorithms[NUM_MH] = all_algorithms[idx];
            NUM_MH++;
        }

        if (strncmp(line, "MAXRUNS", 7) == 0)
        {
            sscanf(line, "MAXRUNS %d", &runData.MAXRUNS);
        }

        else if (strncmp(line, "debug", 5) == 0)
        {
            sscanf(line, "debug %d", &runData.debug);
        }

        else if (strncmp(line, "control", 7) == 0)
        {
            sscanf(line, "control %d", &runData.control);
        }

        else if (strncmp(line, "strategy", 8) == 0)
        {
            sscanf(line, "strategy %d", &runData.strategy);
        }

        else if (strncmp(line, "restart", 7) == 0)
        {
            sscanf(line, "restart %f", &runData.restart);
        }

        else if (strncmp(line, "sizePool", 8) == 0)
        {
            sscanf(line, "sizePool %d", &runData.sizePool);
        }
    }
    fclose(fileConf);

    TProblemData data;
    ReadData(nameInstance, data);

    double foBest = INFINITY,
           f1Best = INFINITY,
           f2Best = INFINITY,
           foAverage = 0.0;

    float timeBest = 0.0,
          timeTotal = 0.0;

    std::vector<double> ofvs;
    ofvs.clear();

    TSol bestSolution;
    bestSolution.ofv = std::numeric_limits<double>::infinity();

    printf("\n\nInstance: %s \nRun: ", nameInstance);
    for (int run = 0; run < runData.MAXRUNS; run++)
    {

        int RSEED = 0;

        RSEED = std::chrono::steady_clock::now().time_since_epoch().count();

        if (!runData.debug)
            rng.seed(RSEED);

        if (runData.debug)
            rng.seed(1234);

        printf("%d ", run + 1);

        TSol bestSolutionRun;
        double f1, f2;
        bestSolution.best_time = 0.0;

        double start_time, end_time;
        start_time = get_time_in_seconds();
        end_time = get_time_in_seconds();

        stop_execution.store(false);

        if (NUM_MH > 0)
        {

            pool.clear();
            pool.resize(runData.sizePool);
            CreatePoolSolutions(data, runData.sizePool);

            bestSolutionRun = pool[0];
            f1 = 0;
            f2 = 0;

            omp_set_num_threads(NUM_MH);
#pragma omp parallel private(rng) shared(pool, stop_execution)
            {
                while (end_time - start_time < runData.MAXTIME)
                {

                    stop_execution.store(false);

#pragma omp for
                    for (int i = 0; i < NUM_MH; ++i)
                    {

#pragma omp cancellation point for

                        void (*function_mh)(const TRunData &, const TProblemData &) = functions_MH[i];

                        if (runData.debug)
                        {
                            printf("\nThread %d executing MH_%d %s [%.2lf].", omp_get_thread_num(), i, algorithms[i], end_time - start_time);
                        }

                        function_mh(runData, data);

                        stop_execution.store(true);
#pragma omp cancel for
                    }

                    end_time = get_time_in_seconds();

#pragma omp critical
                    {
                        if (pool[0].ofv < bestSolutionRun.ofv)
                        {
                            bestSolutionRun = pool[0];
                            f1 = bestSolution.f1;
                            f2 = bestSolution.f2;
                        }
                    }

                    if (end_time - start_time < runData.MAXTIME)
                    {
#pragma omp single
                        {
                            CreatePoolSolutions(data, runData.sizePool);
                        }
                    }
                }
            }

            stop_execution.store(false);
        }

        else
        {
            printf("\n\nThis solver has not been implemented.\n");
            exit(1);
        }

        if (bestSolutionRun.ofv < bestSolution.ofv)
            bestSolution = bestSolutionRun;

        if (bestSolutionRun.ofv < foBest)
        {
            foBest = bestSolutionRun.ofv;
            f1Best = bestSolutionRun.f1;
            f2Best = bestSolutionRun.f2;
        }

        foAverage += bestSolutionRun.ofv;

        ofvs.push_back(bestSolutionRun.ofv);

        timeBest += bestSolutionRun.best_time - start_time;
        timeTotal += end_time - start_time;
    }

    foAverage = foAverage / runData.MAXRUNS;
    timeBest = timeBest / runData.MAXRUNS;
    timeTotal = timeTotal / runData.MAXRUNS;

    if (!runData.debug)
    {
        WriteSolution(algorithms, NUM_MH, bestSolution, timeBest, timeTotal, nameInstance, data);
        WriteResults(algorithms, NUM_MH, bestSolution.ofv, bestSolution.f1, bestSolution.f2,
                     foAverage, ofvs, timeBest, timeTotal, nameInstance, data, bestSolution);
    }
    else
    {
        WriteSolutionScreen(algorithms, NUM_MH, bestSolution, timeBest, timeTotal, nameInstance, data, pool);
    }

    FreeMemoryProblem(data);

    return 0;
}
