#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "Data.h"

const double ESTIMATIVA_E = 4.0;

void WriteSolutionScreen(const char *algorithms[], int numMH, TSol s,
                          float timeBest, float timeTotal, char instance[],
                          const TProblemData &data, std::vector <TSol> pool)
{
    printf("\n\n\nRKO: ");
    for (int i=0; i<numMH; i++)
        printf("%s | ", algorithms[i]);
    printf("\nBest MH: %s \nInstance: %s \nsol: ", s.nameMH, instance);
    for (int i=0; i<data.n; i++)
        printf("%.3lf ", s.rk[i]);

    printf("\nofv: %.5lf", s.ofv);
    printf("\nf1: %.5lf", s.f1);
    printf("\nf2: %.5lf", s.f2);
    printf("\nTotal time: %.3f",timeTotal);
    printf("\nBest time: %.3f\n\n",timeBest);

    printf("\nSolution Pool:\n");
    for (int i = 0; i< (int)pool.size(); i++)
        printf("%.5lf [%s]\n", pool[i].ofv, pool[i].nameMH);

    printf("\n===========================================================\n");
    printf(" RESUMO DA MISSÃO");
    printf("\n===========================================================\n");

    TSol sPrint = s;
    Decoder(sPrint, data);

    double total_victims_found_global = 0.0;
    int total_visits_global = 0;
    double total_energy_consumed_global = 0.0;

    for (int k = 0; k < (int)sPrint.routes.size(); k++)
    {
        if (sPrint.routes[k].size() > 1)
        {
            int start_node_idx = sPrint.routes[k][0];
            double current_battery = data.node[start_node_idx].e;
            double initial_battery = current_battery;
            int visits_count = 0;
            double drone_victims_found = 0.0;
            double energy_consumed = 0.0;

            for (size_t i = 0; i < sPrint.routes[k].size() - 1; i++)
            {
                int u_idx = sPrint.routes[k][i];
                int v_idx = sPrint.routes[k][i+1];

                double dist = data.dist[u_idx][v_idx];
                current_battery -= dist;
                energy_consumed += dist;
                visits_count++;

                double prob = data.node[v_idx].p;
                drone_victims_found += (prob * ESTIMATIVA_E);
            }

            total_victims_found_global += drone_victims_found;
            total_visits_global += visits_count;
            total_energy_consumed_global += energy_consumed;

            printf("DRONE %d: %3d visitas | Vitimas Est.: %6.2lf | Bat. Restante: %6.2lf (%.1f%% usada)\n",
                   k,
                   visits_count,
                   drone_victims_found,
                   current_battery,
                   (energy_consumed/initial_battery)*100.0);
        }
    }
    printf("-----------------------------------------------------------\n");
    printf(" TOTAL GLOBAL VISITAS: %d\n", total_visits_global);
    printf(" TOTAL GLOBAL VÍTIMAS: %.2lf\n", total_victims_found_global);
    printf(" TOTAL ENERGIA GASTA:  %.2lf\n", total_energy_consumed_global);
    printf("===========================================================\n");
     for (int k = 0; k < (int)sPrint.routes.size(); k++)
    {

        if (sPrint.routes[k].size() > 1)
        {
            for (size_t i = 0; i < sPrint.routes[k].size() - 1; i++)
            {
                int u_idx = sPrint.routes[k][i];
                int v_idx = sPrint.routes[k][i+1];

                printf("Drone %d viaja do no %d para o no %d\n",
                       k,
                       data.node[u_idx].id,
                       data.node[v_idx].id);
            }
        }
    }
}

void WriteSolution(const char *algorithms[], int numMH, TSol s,
                    float timeBest, float timeTotal, char instance[],
                    const TProblemData &data)
{
    char name[256]="../Results/Solutions_RKO";
    strcat(name,".txt");

    FILE *solFile = fopen(name,"a");
    if (!solFile) { printf("\n\nFile not found %s!!!",name); getchar(); exit(1); }

    fprintf(solFile,"\n\nInstance: %s", instance);
    fprintf(solFile,"\nRKO: ");
    for (int i=0; i<numMH; i++) fprintf(solFile,"%s | ", algorithms[i]);

    fprintf(solFile,"\nSol: ");
    for (int i=0; i<data.n; i++) fprintf(solFile,"%.3lf ", s.rk[i]);

    fprintf(solFile,"\nofv: %lf", s.ofv);
    fprintf(solFile,"\nf1: %.5lf", s.f1);
    fprintf(solFile,"\nf2: %.5lf", s.f2);
    fprintf(solFile,"\nBest time: %.3f",timeBest);
    fprintf(solFile,"\nTotal time:%.3f \n",timeTotal);

    TSol sFile = s;
    Decoder(sFile, data);

    double total_victims_found_global = 0.0;
    int total_visits_global = 0;
    double total_energy_consumed_global = 0.0;

    fprintf(solFile, "\n--- RESUMO DA MISSÃO ---\n");

    for (int k = 0; k < (int)sFile.routes.size(); k++)
    {
        if (sFile.routes[k].size() > 1)
        {
            int start_node_idx = sFile.routes[k][0];
            double current_battery = data.node[start_node_idx].e;
            double initial_battery = current_battery;
            int visits_count = 0;
            double drone_victims_found = 0.0;
            double energy_consumed = 0.0;

            for (size_t i = 0; i < sFile.routes[k].size() - 1; i++)
            {
                int u_idx = sFile.routes[k][i];
                int v_idx = sFile.routes[k][i+1];

                double dist = data.dist[u_idx][v_idx];
                current_battery -= dist;
                energy_consumed += dist;
                visits_count++;

                double prob = data.node[v_idx].p;
                drone_victims_found += (prob * ESTIMATIVA_E);
            }

            total_victims_found_global += drone_victims_found;
            total_visits_global += visits_count;
            total_energy_consumed_global += energy_consumed;

            fprintf(solFile, "DRONE %d: %d visitas, Vitimas: %.2lf, Bat. Final: %.2lf\n",
                   k, visits_count, drone_victims_found, current_battery);
        }
    }
    fprintf(solFile, "--- TOTAIS: Visitas: %d | Vitimas: %.2lf | Energia: %.2lf ---\n",
            total_visits_global, total_victims_found_global, total_energy_consumed_global);

    fclose(solFile);
}

void WriteResults(const char *algorithms[], int numMH, double ofv, double f1, double f2,
                  double ofvAverage, std::vector <double> ofvs, float timeBest,
                  float timeTotal, char instance[], const TProblemData &data, TSol s)
{
    char name[256]="../Results/Results_RKO";
    strcat(name,".csv");

    FILE *File = fopen(name,"a");
    if (!File) { printf("\n\nFile not found %s!!!",name); exit(1); }

    TSol sCSV = s;
    Decoder(sCSV, data);

    double total_victims_est = 0.0;
    int total_visits = 0;
    double total_energy_consumed = 0.0;

    for (const auto& route : sCSV.routes)
    {
        if (route.size() > 1)
        {
            for (size_t i = 0; i < route.size() - 1; i++)
            {
                int u = route[i];
                int v = route[i+1];

                double dist = data.dist[u][v];
                total_energy_consumed += dist;
                total_visits++;

                double prob = data.node[v].p;
                total_victims_est += (prob * ESTIMATIVA_E);
            }
        }
    }

    fprintf(File,"\n%s\t", instance);
    for (int i=0; i<numMH; i++) fprintf(File,"%s | ", algorithms[i]);

    fprintf(File,"\t%lf", ofv);
    fprintf(File,"\t%lf", f1);
    fprintf(File,"\t%lf", f2);
    fprintf(File,"\t%lf", ofvAverage);
    fprintf(File,"\t%.3f", timeBest);
    fprintf(File,"\t%.3f", timeTotal);

    fprintf(File,"\t%d", total_visits);
    fprintf(File,"\t%.4lf", total_victims_est);
    fprintf(File,"\t%.4lf", total_energy_consumed);

    fprintf(File,"\t%d", (int)ofvs.size());
    for (unsigned int i=0; i<ofvs.size(); i++){
        fprintf(File,"\t%lf", ofvs[i]);
    }

    fclose(File);
}
