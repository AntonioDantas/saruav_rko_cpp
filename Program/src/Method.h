#ifndef _Method_H
#define _Method_H

bool sortByFitness(const TSol &lhs, const TSol &rhs) { return lhs.ofv < rhs.ofv; }

double randomico(double min, double max)
{
    double val;
#pragma omp critical(rng_lock)
    {
        val = std::uniform_real_distribution<double>(min, max)(rng);
    }
    return val;
}

int irandomico(int min, int max)
{
    return (int)randomico(0, max - min + 1) + min;
}

double get_time_in_seconds()
{
#if defined(_WIN32) || defined(_WIN64)
    LARGE_INTEGER frequency, timeCur;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&timeCur);
    return static_cast<double>(timeCur.QuadPart) / frequency.QuadPart;
#else
    struct timespec timeCur;
    clock_gettime(CLOCK_MONOTONIC, &timeCur);
    return timeCur.tv_sec + timeCur.tv_nsec / 1e9;
#endif
}

void CreateInitialSolutions(TSol &s, const int n)
{
    s.rk.resize(n);

    for (int j = 0; j < n; j++)
    {
        s.rk[j] = randomico(0, 1);
    }
}

void CreatePoolSolutions(const TProblemData &data, const int sizePool)
{
#pragma omp critical
    {
        for (int i = 0; i < sizePool; i++)
        {
            CreateInitialSolutions(pool[i], data.n);
            pool[i].ofv = Decoder(pool[i], data);
            pool[i].best_time = get_time_in_seconds();
        }

        for (int i = 0; i < sizePool - 1; i++)
        {
            int minIndex = i;

            for (int j = i + 1; j < sizePool; j++)
            {
                if (pool[j].ofv < pool[minIndex].ofv)
                {
                    minIndex = j;
                }
            }

            if (minIndex != i)
            {
                TSol aux = pool[i];
                pool[i] = pool[minIndex];
                pool[minIndex] = aux;
            }
        }

        int clone = 0;
        for (int i = sizePool - 1; i > 0; i--)
        {
            if (pool[i].ofv == pool[i - 1].ofv)
            {
                for (int j = 0; j < 0.2 * data.n; j++)
                {
                    int pos = irandomico(0, data.n - 1);
                    pool[i].rk[pos] = randomico(0, 1);
                }
                pool[i].ofv = Decoder(pool[i], data);
                clone = 1;
            }
        }

        if (clone)
        {
            for (int i = 0; i < sizePool - 1; i++)
            {
                int minIndex = i;

                for (int j = i + 1; j < sizePool; j++)
                {
                    if (pool[j].ofv < pool[minIndex].ofv)
                    {
                        minIndex = j;
                    }
                }

                if (minIndex != i)
                {
                    TSol aux = pool[i];
                    pool[i] = pool[minIndex];
                    pool[minIndex] = aux;
                }
            }
        }
    }
}

void UpdatePoolSolutions(TSol s, const char *mh, const int debug)
{
#pragma omp critical
    {

        bool exists = false;
        for (int i = 0; i < (int)pool.size(); i++)
        {
            if (pool[i].ofv == s.ofv)
            {
                exists = true;
                break;
            }
        }

        if (s.ofv < pool[0].ofv && debug)
        {

            int thread_id = omp_get_thread_num();
            printf("\nBest solution: %.10lf (Thread: %d - MH: %s)", s.ofv, thread_id, mh);
        }

        if (!exists)
        {

            s.best_time = get_time_in_seconds();

            strcpy(s.nameMH, mh);

            int i;
            for (i = (int)pool.size() - 1; i > 0 && pool[i - 1].ofv > s.ofv; i--)
            {
                pool[i] = pool[i - 1];
            }

            pool[i] = s;
        }
    }
}

void ShakeSolution(TSol &s, float betaMin, float betaMax, const int n)
{
    int shaking_type = 0;
    int intensity = (int)(n * randomico(betaMin, betaMax)) + 1;
    if (intensity < 1)
        intensity = 1;
    for (int k = 0; k < intensity; k++)
    {
        shaking_type = irandomico(1, 4);
        int i = irandomico(0, n - 1);

        if (shaking_type == 1)
        {

            s.rk[i] = randomico(0, 1);
        }
        else if (shaking_type == 2)
        {

            if (s.rk[i] > 0.0001)
                s.rk[i] = 1.0 - s.rk[i];
            else
                s.rk[i] = 0.9999;
        }
        else if (shaking_type == 3)
        {

            int j = irandomico(0, n - 1);
            double temp = s.rk[i];
            s.rk[i] = s.rk[j];
            s.rk[j] = temp;
        }
        i = irandomico(0, n - 2);
        if (shaking_type == 4)
        {

            double temp = s.rk[i];
            s.rk[i] = s.rk[i + 1];
            s.rk[i + 1] = temp;
        }
    }
}

TSol Blending(TSol &s1, TSol &s2, double factor, const int n)
{
    TSol s;

    s.rk.clear();
    s.rk.resize(n);

    for (int j = 0; j < n; j++)
    {

        if (randomico(0, 1) < 0.02)
        {
            s.rk[j] = randomico(0, 1);
        }

        else
        {
            if (randomico(0, 1) < 0.5)
            {
                s.rk[j] = s1.rk[j];
            }
            else
            {
                if (factor == -1)
                {
                    s.rk[j] = std::clamp(1.0 - s2.rk[j], 0.0, 0.9999999);
                }
                else
                {
                    s.rk[j] = s2.rk[j];
                }
            }
        }
    }
    return s;
}

void NelderMeadSearch(TSol &x1, const TProblemData &data)
{
    int improved = 0;
    int improvedX1 = 0;
    TSol x1Origem = x1;

    int k1, k2;
    TSol x2, x3;

#pragma omp critical
    {
        do
        {
            k1 = irandomico(0, pool.size() - 1);
            k2 = irandomico(0, pool.size() - 1);
        } while (k1 == k2);

        x2 = pool[k1];
        x3 = pool[k2];
    }

    TSol x_r;
    TSol x_e;
    TSol x_c;
    TSol x0;

    TSol xBest = x1;

    int iter_count = 0;
    int eval_count = 0;

    if (x1.ofv > x2.ofv)
    {
        TSol temp = x1;
        x1 = x2;
        x2 = temp;
    }

    if (x1.ofv > x3.ofv)
    {
        TSol temp = x1;
        x1 = x3;
        x3 = temp;
    }

    if (x2.ofv > x3.ofv)
    {
        TSol temp = x2;
        x2 = x3;
        x3 = temp;
    }

    x0 = Blending(x1, x2, 1, data.n);
    x0.ofv = Decoder(x0, data);
    if (x0.ofv < xBest.ofv)
    {
        xBest = x0;
        improved = 1;
    }

    iter_count++;

    int maxIter = data.n * exp(-2);
    while (iter_count <= maxIter)
    {
        int shrink = 0;

        x_r = Blending(x0, x3, -1, data.n);
        x_r.ofv = Decoder(x_r, data);
        if (x_r.ofv < xBest.ofv)
        {
            xBest = x_r;
            improved = 1;
            improvedX1 = 1;
        }
        eval_count++;

        if (x_r.ofv < x1.ofv)
        {

            x_e = Blending(x_r, x0, -1, data.n);
            x_e.ofv = Decoder(x_e, data);
            if (x_e.ofv < xBest.ofv)
            {
                xBest = x_e;
                improved = 1;
                improvedX1 = 1;
            }
            eval_count++;

            if (x_e.ofv < x_r.ofv)
            {

                x3 = x_e;
            }
            else
            {

                x3 = x_r;
            }
        }

        else
        {

            if (x_r.ofv < x2.ofv)
            {

                x3 = x_r;
            }
            else
            {

                if (x_r.ofv < x3.ofv)
                {

                    x_c = Blending(x_r, x0, 1, data.n);
                    x_c.ofv = Decoder(x_c, data);
                    if (x_c.ofv < xBest.ofv)
                    {
                        xBest = x_c;
                        improved = 1;
                        improvedX1 = 1;
                    }
                    eval_count++;

                    if (x_c.ofv < x_r.ofv)
                    {

                        x3 = x_c;
                    }
                    else
                    {

                        shrink = 1;
                    }
                }
                else
                {

                    x_c = Blending(x0, x3, 1, data.n);
                    x_c.ofv = Decoder(x_c, data);
                    if (x_c.ofv < xBest.ofv)
                    {
                        xBest = x_c;
                        improved = 1;
                        improvedX1 = 1;
                    }
                    eval_count++;

                    if (x_c.ofv < x3.ofv)
                    {

                        x3 = x_c;
                    }
                    else
                    {

                        shrink = 1;
                    }
                }
            }
        }
        if (shrink)
        {
            x2 = Blending(x1, x2, 1, data.n);
            x2.ofv = Decoder(x2, data);
            if (x2.ofv < xBest.ofv)
            {
                xBest = x2;
                improved = 1;
                improvedX1 = 1;
            }

            eval_count++;

            x3 = Blending(x1, x3, 1, data.n);
            x3.ofv = Decoder(x3, data);
            if (x3.ofv < xBest.ofv)
            {
                xBest = x3;
                improved = 1;
                improvedX1 = 1;
            }
            eval_count++;
        }

        if (x1.ofv > x2.ofv)
        {
            TSol temp = x1;
            x1 = x2;
            x2 = temp;
        }

        if (x1.ofv > x3.ofv)
        {
            TSol temp = x1;
            x1 = x3;
            x3 = temp;
        }

        if (x2.ofv > x3.ofv)
        {
            TSol temp = x2;
            x2 = x3;
            x3 = temp;
        }

        x0 = Blending(x1, x2, 1, data.n);
        x0.ofv = Decoder(x0, data);
        if (x0.ofv < xBest.ofv)
        {
            xBest = x0;
            improved = 1;
            improvedX1 = 1;
        }

        if (improved == 1)
        {
            improved = 0;
            iter_count = 0;
        }
        else
            iter_count++;

        if (stop_execution.load())
            return;
    }

    if (improvedX1 == 1)
    {

        x1 = xBest;
    }
    else
    {

        x1 = x1Origem;
    }
}

void SwapLS(TSol &s, const TProblemData &data, const int &strategy, std::vector<int> &RKorder)
{

    std::shuffle(RKorder.begin(), RKorder.end(), rng);

    float rate = 1.0;

    if (strategy == 1)
    {
        TSol sBest = s;
        for (int i = 0; i < (data.n - 1) * rate; i++)
        {
            for (int j = i + 1; j < data.n * rate; j++)
            {

                double temp = s.rk[RKorder[i]];
                s.rk[RKorder[i]] = s.rk[RKorder[j]];
                s.rk[RKorder[j]] = temp;

                s.ofv = Decoder(s, data);

                if (s.ofv < sBest.ofv)
                {
                    sBest = s;

                }

                else
                {
                    s = sBest;
                }
                if (stop_execution.load())
                    return;
            }
        }
    }

    if (strategy == 2)
    {
        TSol sBest = s;
        TSol sCurrent = s;
        int improved = 0;
        for (int i = 0; i < (data.n - 1) * rate; i++)
        {
            for (int j = i + 1; j < data.n; j++)
            {

                double temp = s.rk[RKorder[i]];
                s.rk[RKorder[i]] = s.rk[RKorder[j]];
                s.rk[RKorder[j]] = temp;

                s.ofv = Decoder(s, data);
                if (s.ofv < sBest.ofv)
                {
                    sBest = s;
                    improved = 1;
                }

                s = sCurrent;

                if (stop_execution.load())
                    return;
            }
        }
        if (improved == 1)
            s = sBest;
        else
            s = sCurrent;
    }
}

void InvertLS(TSol &s, const TProblemData &data, const int &strategy, std::vector<int> &RKorder)
{

    std::shuffle(RKorder.begin(), RKorder.end(), rng);

    float rate = 1.0;

    if (strategy == 1)
    {
        TSol sBest = s;
        for (int i = 0; i < data.n * rate; i++)
        {

            if (s.rk[RKorder[i]] > 0.00001)
                s.rk[RKorder[i]] = 1.0 - s.rk[RKorder[i]];
            else
                s.rk[RKorder[i]] = 0.99999;

            s.ofv = Decoder(s, data);
            if (s.ofv < sBest.ofv)
            {
                sBest = s;

            }
            else
            {
                s = sBest;
            }
            if (stop_execution.load())
                return;
        }
    }

    else
    {
        TSol sBest = s;
        TSol sCurrent = s;
        int improved = 0;
        for (int i = 0; i < data.n * rate; i++)
        {

            if (s.rk[RKorder[i]] > 0.00001)
                s.rk[RKorder[i]] = 1.0 - s.rk[RKorder[i]];
            else
                s.rk[RKorder[i]] = 0.99999;

            s.ofv = Decoder(s, data);
            if (s.ofv < sBest.ofv)
            {
                sBest = s;
                improved = 1;
            }

            s = sCurrent;

            if (stop_execution.load())
                return;
        }

        if (improved == 1)
            s = sBest;
        else
            s = sCurrent;
    }
}

void FareyLS(TSol &s, const TProblemData &data, const int &strategy, std::vector<int> &RKorder)
{

    std::shuffle(RKorder.begin(), RKorder.end(), rng);

    std::vector<double> F = {0.00, 0.142857, 0.166667, 0.20, 0.25, 0.285714, 0.333333, 0.40, 0.428571, 0.50,
                             0.571429, 0.60, 0.666667, 0.714286, 0.75, 0.80, 0.833333, 0.857143, 1.0};

    float rate = 1.0;

    if (strategy == 1)
    {
        TSol sBest = s;
        for (int i = 0; i < data.n * rate; i++)
        {
            for (int j = 0; j < (int)F.size() - 1; j++)
            {

                s.rk[RKorder[i]] = randomico(F[j], F[j + 1]);

                s.ofv = Decoder(s, data);

                if (s.ofv < sBest.ofv)
                {
                    sBest = s;

                }
                else
                {
                    s = sBest;
                }
                if (stop_execution.load())
                    return;
            }
        }
    }

    else
    {
        TSol sBest = s;
        TSol sCurrent = s;
        int improved = 0;

        for (int i = 0; i < data.n * rate; i++)
        {
            for (int j = 0; j < (int)F.size() - 1; j++)
            {

                s.rk[RKorder[i]] = randomico(F[j], F[j + 1]);

                s.ofv = Decoder(s, data);
                if (s.ofv < sBest.ofv)
                {
                    sBest = s;
                    improved = 1;
                }

                s = sCurrent;
                if (stop_execution.load())
                    return;
            }
        }
        if (improved == 1)
            s = sBest;
        else
            s = sCurrent;
    }
}

void RVND(TSol &s, const TProblemData &data, const int &strategy, std::vector<int> &RKorder)
{

    int numLS = 4;

    int k = 1;

    std::vector<int> NSL;
    std::vector<int> NSLAux;

    for (int i = 1; i <= numLS; i++)
    {
        NSL.push_back(i);
        NSLAux.push_back(i);
    }

    int numIter = 0;

    while (!NSL.empty())
    {
        if (stop_execution.load())
            return;

        double foCurrent = s.ofv;
        numIter++;

        int pos = rand() % NSL.size();
        k = NSL[pos];

        switch (k)
        {
        case 1:
            SwapLS(s, data, strategy, RKorder);
            break;

        case 2:
            InvertLS(s, data, strategy, RKorder);
            break;

        case 3:
            NelderMeadSearch(s, data);
            break;

        case 4:
            FareyLS(s, data, strategy, RKorder);
            break;

        default:
            break;
        }

        if (s.ofv < foCurrent)
        {

            NSL.clear();
            NSL = NSLAux;
        }
        else
        {

            NSL.erase(NSL.begin() + pos);
        }
    }
}

void readParameters(const char *method, int control,
                    std::vector<std::vector<double>> &parameters, int numPar)
{
#pragma omp critical
    {
        char paramFile[256];
        if (control == 0)
        {
            strncpy(paramFile, "config/ParametersOffline.txt", 255);
        }
        else
        {
            strncpy(paramFile, "config/ParametersOnline.txt", 255);
        }

        FILE *file = fopen(paramFile, "r");
        if (file == NULL)
        {
            printf("Error in open file %s.\n", paramFile);
            getchar();
        }

        char line[100];

        while (fgets(line, sizeof(line), file) != NULL)
        {
            line[strcspn(line, "\n")] = '\0';

            if (strcmp(line, method) == 0)
            {

                for (int i = 0; i < numPar; i++)
                {
                    double aux = 0;
                    if (fgets(line, sizeof(line), file) != NULL)
                    {
                        char *token = strtok(line, "{},= ");
                        while (token != NULL)
                        {
                            if (sscanf(token, "%lf", &aux) == 1)
                            {
                                parameters[i].push_back(aux);
                            }
                            token = strtok(NULL, "{},= ");
                        }
                    }
                }
            }
        }
        fclose(file);
    }
}

#endif
