#ifndef _PSO_H
#define _PSO_H

static void UpdateParticleSize(std::vector <TSol> &X, std::vector <TSol> &Pbest, std::vector <std::vector <float> > &V, TSol Gbest, int Psize,
                              float c1, float c2, float w, const TProblemData &data)
{

    int oldPsize = X.size();

    if (oldPsize > Psize){
        X.resize(Psize);
        Pbest.resize(Psize);
        V.resize(Psize);
    }

    else if (oldPsize < Psize){
        X.resize(Psize);
        Pbest.resize(Psize);
        V.resize(Psize, std::vector<float>(data.n));

        for (int i=oldPsize; i<Psize; i++)
        {

            CreateInitialSolutions(X[i], data.n);

            Pbest[i] = X[i];

            for (int j=0; j<data.n; j++)
                V[i][j] = randomico(0,1);

            for (int j=0; j<data.n; j++)
            {
                float r1 = randomico(0,1);
                float r2 = randomico(0,1);

                V[i][j] = w * (V[i][j] + (c1 * r1 * (Pbest[i].rk[j] - X[i].rk[j])) +
                                         (c2 * r2 * (Gbest.rk[j] - X[i].rk[j])));

                double oldrk = X[i].rk[j];
                X[i].rk[j] = X[i].rk[j] + V[i][j];

                if (X[i].rk[j] < 0.0 || X[i].rk[j] >= 1.0) {
                    X[i].rk[j] = oldrk;
                    V[i][j] = 0;
                }
            }

            X[i].ofv = Decoder(X[i], data);

            Pbest[i] = X[i];
        }
    }
}

void PSO(const TRunData &runData, const TProblemData &data)
{
    const char* method = "PSO";
    int Psize = 0;
    float c1 = 0.0;
    float c2 = 0.0;
    float w = 0.0;

    std::vector <TSol> X;
    std::vector <TSol> Pbest;
    std::vector <std::vector <float> > V;

    TSol Gbest;
    double bestOFcurrent = 0;

    int numGenerations = 0;
    float currentTime = 0;
    int bestGeneration = 0;
    int improv = 0;

    double start_timeMH = get_time_in_seconds();
    double end_timeMH = get_time_in_seconds();

    std::vector<int> RKorder;
    RKorder.resize(data.n);
    std::iota(RKorder.begin(), RKorder.end(), 0);

    std::vector<TState> S;
    int numPar = 0;
    int numStates = 0;
    int iCurr = 0;
    double epsilon=0;
    double lf=0;
    double df=0;
    double R=0;
    std::vector <std::vector <TQ> > Q;
    std::vector<int> ai;
    float epsilon_max = 1.0;
    float epsilon_min = 0.1;
    int Ti = 1;
    int restartEpsilon = 1;
    int st = 0;
    int at = 0;

    numPar = 4;
    std::vector<std::vector<double>> parameters;
    parameters.resize(numPar);

    readParameters(method, runData.control, parameters, numPar);

    if (runData.control == 0){

        Psize = parameters[0][0];
        c1    = parameters[1][0];
        c2    = parameters[2][0];
        w     = parameters[3][0];
    }

    else{

        if (runData.control == 1){

            CreateStates(parameters, numStates, numPar, S);

            restartEpsilon = 1;

            epsilon_max = 1.0;

            iCurr = irandomico(0,numStates-1);

            Psize = (int)S[iCurr].par[0];
            c1 = S[iCurr].par[1];
            c2 = S[iCurr].par[2];
            w  = S[iCurr].par[3];
        }
    }

    X.clear();
    Pbest.clear();
    V.clear();

    X.resize(Psize);
    Pbest.resize(Psize);
    V.resize(Psize, std::vector<float>(data.n));

    Gbest.ofv = INFINITY;
    for (int i=0; i<Psize; i++)
    {

        CreateInitialSolutions(X[i], data.n);

        X[i].ofv = Decoder(X[i], data);

        for (int j=0; j<data.n; j++)
            V[i][j] = randomico(0,1);

        if (X[i].ofv < Gbest.ofv)
            Gbest = X[i];

        Pbest[i] = X[i];
    }

    while (currentTime < runData.MAXTIME*runData.restart)
    {

        numGenerations++;

        if (runData.control == 1){

            SetQLParameter(currentTime, Ti, restartEpsilon, epsilon_max, epsilon_min, epsilon, lf, df, runData.MAXTIME*runData.restart);

            at = ChooseAction(S, st, epsilon);

            iCurr = S[st].Ai[at];

            Psize = (int)S[iCurr].par[0];
            c1 = S[iCurr].par[1];
            c2 = S[iCurr].par[2];
            w  = S[iCurr].par[3];

            UpdateParticleSize(X, Pbest, V, Gbest, Psize, c1, c2, w, data);
        }

        bestOFcurrent = INFINITY;
        double media = 0;
        for (int i=0; i<Psize; i++)
        {
            if (stop_execution.load()) return;

            double probUpdate = 1.0;

            for (int j=0; j<data.n; j++)
            {
                float r1 = randomico(0,1);
                float r2 = randomico(0,1);

                V[i][j] = w * (V[i][j] + (c1 * r1 * (Pbest[i].rk[j] - X[i].rk[j])) +
                                         (c2 * r2 * (Gbest.rk[j] - X[i].rk[j])));

                if (randomico(0,1) < probUpdate){

                    double oldrk = X[i].rk[j];
                    X[i].rk[j] = X[i].rk[j] + V[i][j];

                    if (X[i].rk[j] < 0.0 || X[i].rk[j] >= 1.0) {
                        X[i].rk[j] = oldrk;
                        V[i][j] = 0;
                    }
                }
            }

            X[i].ofv = Decoder(X[i], data);

            if (X[i].ofv < bestOFcurrent){
                bestOFcurrent = X[i].ofv;
            }

            if (X[i].ofv < Pbest[i].ofv){
                Pbest[i] = X[i];
            }

            if (X[i].ofv < Gbest.ofv){
                Gbest = X[i];
                bestGeneration = numGenerations;
                improv = 1;
            }

            media += X[i].ofv;
        }

        double oldGbest = Gbest.ofv;
        int chosen = irandomico(0,Psize-1);
        NelderMeadSearch(Pbest[chosen], data);

        if (Pbest[chosen].ofv < Gbest.ofv){
            Gbest = Pbest[chosen];
            bestGeneration = numGenerations;
            improv = 1;
        }

        if (bestGeneration == numGenerations || Gbest.ofv < oldGbest){

            UpdatePoolSolutions(Gbest, method, runData.debug);
        }

        if (runData.control == 1){

            if (improv){
                R = 1;
                improv = 0;
            }
            else{
                R = (Gbest.ofv - bestOFcurrent)/bestOFcurrent;
            }

            int st_1 = S[st].Ai[at];

            S[st].Qa[at] = S[st].Qa[at] + lf*(R + df*S[st_1].maxQ - S[st].Qa[at]);

            if (S[st].Qa[at] > S[st].maxQ)
            {
                S[st].maxQ = S[st].Qa[at];
                S[st].maxA = at;
            }

            st = st_1;
        }

        end_timeMH = get_time_in_seconds();
        currentTime = end_timeMH - start_timeMH;
    }

    X.clear();
    V.clear();
    Pbest.clear();

}

#endif
