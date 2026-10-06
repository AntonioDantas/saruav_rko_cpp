#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Data.h"
#include <ranges>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <limits>

struct TNode
{
    int id;
    double p;
    double e;
};

struct TProblemData
{
    int n;

    int nodeCount;
    int veichesCount;
    int vMin;
    int vMax;
    double f1Max;
    double f2Max;
    double total_probability;
    double max_potential_score;

    std::vector<TNode> node;
    std::vector<std::vector<double>> dist;
    std::vector<std::vector<double>> distNorm;
};

const double PENALTY_FACTOR = 99999.0;

double find_max_matrix(const std::vector<std::vector<double>> &m, size_t nRowsCols)
{
    double maxV = -std::numeric_limits<double>::infinity();
    size_t n = nRowsCols;
    for (size_t i = 0; i < n && i < m.size(); ++i)
    {
        for (size_t j = 0; j < n && j < m[i].size(); ++j)
        {
            double val = m[i][j];
            if (!std::isfinite(val))
                continue;
            if (std::fabs(val) <= 0.0)
                continue;
            if (val > maxV)
                maxV = val;
        }
    }
    return maxV;
}

void normalize_matrix(std::vector<std::vector<double>> &m, size_t nRowsCols)
{
    auto maxV = find_max_matrix(m, nRowsCols);

    if (!std::isfinite(maxV) || maxV <= 0.0)
        return;

    size_t n = nRowsCols;
    for (size_t i = 0; i < n && i < m.size(); ++i)
    {
        for (size_t j = 0; j < n && j < m[i].size(); ++j)
        {
            double val = m[i][j];
            if (!std::isfinite(val) || std::fabs(val) <= 0.0)
                continue;
            m[i][j] = val / maxV;
        }
    }
}

double calculate_route_distance(const std::vector<int> &route, const std::vector<std::vector<double>> &dist)
{
    double total_distance = 0.0;
    if (route.size() < 2)
        return 0.0;

    for (size_t i = 0; i < route.size() - 1; ++i)
    {
        int from_node = route[i];
        int to_node = route[i + 1];
        total_distance += dist[from_node][to_node];
    }
    return total_distance;
}

void ReadData(char nameTable[], TProblemData &data)
{
    data.nodeCount = 0;
    data.veichesCount = 0;
    char name[400] = "../Instances/";
    strcat(name, nameTable);

    FILE *arq;
    arq = fopen(name, "r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n", name);
        getchar();
        exit(1);
    }

    char temp[100];
    fscanf(arq, "%s %d", temp, &data.nodeCount);
    fscanf(arq, "%s %d", temp, &data.veichesCount);
    fscanf(arq, "%s %d", temp, &data.vMin);
    fscanf(arq, "%s %d", temp, &data.vMax);
    fscanf(arq, "%s %lf", temp, &data.f1Max);
    fscanf(arq, "%s %lf", temp, &data.f2Max);

    printf("\nInstance: %s", nameTable);
    printf("\nNumber of Nodes: %d", data.nodeCount);
    printf("\nNumber of Vehicles: %d", data.veichesCount);
    printf("\nMinimum Visits: %d", data.vMin);
    printf("\nMaximum Visits: %d", data.vMax);
    printf("\nMax value of f1: %lf", data.f1Max);
    printf("\nMin value of f2: %lf\n", data.f2Max);

    fscanf(arq, "%s", temp);

    data.node.clear();
    TNode nodeTemp;

    while (data.node.size() < (size_t)data.nodeCount)
    {
        fscanf(arq, "%d %lf %lf", &nodeTemp.id, &nodeTemp.p, &nodeTemp.e);
        data.node.push_back(nodeTemp);
    }

    data.dist.clear();
    data.distNorm.clear();
    data.dist.resize(data.nodeCount, std::vector<double>(data.nodeCount));
    data.distNorm.resize(data.nodeCount, std::vector<double>(data.nodeCount));

    fscanf(arq, "%s", temp);

    for (int i = 0; i < data.nodeCount; i++)
    {
        for (int j = 0; j < data.nodeCount; j++)
        {
            if (fscanf(arq, "%lf", &data.dist[i][j]) != 1)
            {
                printf("Erro de leitura no elemento [%d][%d]\n", i, j);
            }
            data.distNorm[i][j] = data.dist[i][j];
        }
    }
    fclose(arq);

    data.total_probability = 0.0;
    data.max_potential_score = 0.0;
    for (const auto& node : data.node) {
        if (node.e <= 0) {
            data.total_probability += node.p;
            data.max_potential_score += node.p * (double)data.vMax;
        }
    }
    printf("\nTotal Potential Score: %.3lf\n", data.max_potential_score);
    printf("\nTotal Probability of Target Nodes: %.3lf\n", data.total_probability);

    normalize_matrix(data.distNorm, data.nodeCount);

    data.n = data.nodeCount + (data.nodeCount - data.veichesCount) * 2;
}

double Decoder(TSol &s, const TProblemData &data)
{
    if (s.rk.size() < (size_t)data.n)
    {
        printf("\nERROR: Size of rk vector (%zu) is less than required (%d)!\n", s.rk.size(), data.n);
        return PENALTY_FACTOR;
    }

    for (double val : s.rk)
    {
        if (std::isnan(val)) return PENALTY_FACTOR;
    }

    s.routes.clear();
    s.routes.reserve(data.veichesCount);
    s.ofv = 0.0;
    s.f1 = 0.0;
    s.f2 = 0.0;

    struct DecodedItem {
        int id;
        double sortKey;
    };
    std::vector<DecodedItem> sequence;
    sequence.reserve(data.nodeCount * 2);

    for (int i = 0; i < data.nodeCount; ++i)
    {
        if (data.node[i].e > 0)
        {
            sequence.push_back({i, s.rk[i]});

            s.routes.push_back({i});
        }
    }

    int idx_check_visit = data.nodeCount;
    int idx_extra_visit_count = data.nodeCount;

    int target_count = 0;
    for(int i=0; i<data.nodeCount; ++i) if(data.node[i].e == 0) target_count++;

    idx_extra_visit_count += target_count;

    for (int i = 0; i < data.nodeCount; ++i)
    {
        if (data.node[i].e == 0)
        {

            double visit_key = s.rk[idx_check_visit++];
            bool is_visited = (data.vMin + floor(visit_key * (data.vMax - data.vMin + 1.0))) > 0;

            if (is_visited)
            {

                sequence.push_back({i, s.rk[i]});

                double revisit_key = s.rk[idx_extra_visit_count++];
                int max_possible_extras = data.vMax - 1;
                double revisit_factor = (double)data.vMax;
                if(data.vMin >= 0) revisit_factor++;
                int extra_visits = floor(revisit_key * revisit_factor);
                if (extra_visits > max_possible_extras) {
                    extra_visits = max_possible_extras;
                }

                for (int k = 0; k < extra_visits; ++k)
                {

                    double derived_key = std::fmod(revisit_key + (k + 1) * 0.61803398875, 1.0);

                    sequence.push_back({i, derived_key});
                }
            }
            else
            {

                idx_extra_visit_count++;
            }
        }
    }

    std::sort(sequence.begin(), sequence.end(),
              [](const DecodedItem &a, const DecodedItem &b)
              { return a.sortKey < b.sortKey; });

    int current_vehicle_id = -1;
    std::vector<int> pending_buffer;

    for (const auto &item : sequence)
    {
        if (data.node[item.id].e > 0)
        {

            current_vehicle_id = item.id;

            if (!pending_buffer.empty())
            {
                auto it = std::find_if(s.routes.begin(), s.routes.end(),
                                     [&](const std::vector<int> &r) { return r[0] == current_vehicle_id; });
                for (int pending_node : pending_buffer)
                {
                    if (it->back() != pending_node)
                    {
                        it->push_back(pending_node);
                    }
                }
                pending_buffer.clear();
            }
        }
        else
        {
            if (current_vehicle_id != -1)
            {

                auto it = std::find_if(s.routes.begin(), s.routes.end(),
                                     [&](const std::vector<int> &r) { return r[0] == current_vehicle_id; });
                if (it->back() != item.id)
                {
                    it->push_back(item.id);
                }
            }
            else
            {

                if (pending_buffer.empty() || pending_buffer.back() != item.id)
                {
                    pending_buffer.push_back(item.id);
                }
            }
        }
    }

    if (!pending_buffer.empty() && current_vehicle_id != -1)
    {
        auto it = std::find_if(s.routes.begin(), s.routes.end(),
                             [&](const std::vector<int> &r) { return r[0] == current_vehicle_id; });
        if (it != s.routes.end())
        {
            for (int pending_node : pending_buffer)
            {
                if (it->back() != pending_node)
                {
                    it->push_back(pending_node);
                }
            }
        }
    }

    double total_energy_overflow = 0.0;
    double total_distance = 0.0;

    for (const auto &route : s.routes)
    {
        if (route.size() <= 1)
            continue;

        double current_distance = 0.0;
        double current_distance_norm = 0.0;

        for (size_t i = 1; i < route.size(); ++i)
        {
            int from_node = route[i - 1];
            int to_node = route[i];
            current_distance += data.dist[from_node][to_node];
            current_distance_norm += data.distNorm[from_node][to_node];

            s.f2 += data.node[to_node].p;

        }
        s.f1 += current_distance_norm;
        total_distance += current_distance;

        if (current_distance > data.node[route[0]].e)
        {
            total_energy_overflow += (current_distance - data.node[route[0]].e);
        }

        if (data.f2Max == 0 && data.f1Max > 0 && current_distance > data.f1Max)
        {
            total_energy_overflow += (current_distance - data.f1Max);
        }
    }

    double f2_normalized = 0.0;
    if (data.max_potential_score > 0) {
        f2_normalized = s.f2 / data.max_potential_score;
    }

    if (data.f1Max == 0)
    {
        s.ofv = s.f1;

        if (f2_normalized < data.f2Max) {
            s.ofv += PENALTY_FACTOR * (data.f2Max - f2_normalized);
        }
    }
    else
    {
        s.ofv = -1.0 * f2_normalized;
    }

    if (total_energy_overflow > 0)
    {
        s.ofv += PENALTY_FACTOR * total_energy_overflow;
    }

    return s.ofv;
}

void FreeMemoryProblem(TProblemData &data)
{
    data.dist.clear();
    data.distNorm.clear();
    data.node.clear();
}

#endif
