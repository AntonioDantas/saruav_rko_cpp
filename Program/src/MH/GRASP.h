#ifndef _GRASP_H
#define _GRASP_H

static void LineSearch(TSol s, float h, int i, double &bestZ, double &bestF, const TProblemData &data)
{

    bestZ = 0;
    bestF = INFINITY;

    double tau = 0;
    std::vector<double> rk;
    rk.push_back(s.rk[i] + tau * h);
    for (int j=0; j<(int)(1.0/h)+1; j+=2){
        tau++;

        if (((s.rk[i] + tau * h) >= 0) && ((s.rk[i] + tau * h) < 1))
            rk.push_back(s.rk[i] + tau * h);

        if (((s.rk[i] + (-1*tau) * h) >= 0) && ((s.rk[i] + (-1*tau) * h) < 1))
            rk.push_back(s.rk[i] + (-1*tau) * h);
    }

    int q = ceil(log2((int)(1.0/h))) + 1;
    if (q > (int)rk.size())
        q = rk.size();

    std::shuffle(rk.begin(), rk.end(), std::mt19937(std::random_device()()));

    for (int j=0; j<q; j++)
    {
        if (stop_execution.load()) return;

        s.rk[i] = rk[j];
        s.ofv = Decoder(s, data);

        if (s.ofv < bestF){
            bestZ = s.rk[i];
            bestF = s.ofv;
        }
    }

    rk.clear();
}

static void ConstrutiveGreedyRandomized(TSol &s, float h, float alpha, const TProblemData &data)
{
    std::vector<int> UnFixed(data.n);
    std::vector<int> chosenRK;
    std::vector<int> RCL;
    std::vector<double> z(data.n);
    std::vector<double> g(data.n,INFINITY);

    double min, max;

    double betaMin = 0.3,
           betaMax = 0.7;

    std::iota(UnFixed.begin(), UnFixed.end(), 0);

    double intensity = randomico(betaMin, betaMax);

    int kBest = 0;

    for (int j=0; j<data.n*intensity; j++)
    {

        min = INFINITY;
        max = -INFINITY;

        chosenRK.clear();
        int kMax = UnFixed.size() * 0.1;
        if (kMax < 2)
            kMax = UnFixed.size();

        chosenRK = UnFixed;
        std::shuffle(chosenRK.begin(), chosenRK.end(), rng);
        chosenRK.resize(kMax);

        z.clear();
        z.resize(data.n);
        g.clear();
        g.resize(data.n,INFINITY);

        for (int k=0; k<kMax; k++)
        {
            int i = chosenRK[k];

            TSol sAux = s;

            LineSearch(sAux, h, i, z[i], g[i], data);

            if (min > g[i])
            {
                min = g[i];
                kBest = i;
            }
        }

        s.rk[kBest] = z[kBest];
        s.ofv = g[kBest];

        for (int i=0; i<(int)UnFixed.size(); i++)
        {
            if (UnFixed[i] == kBest)
            {
                UnFixed.erase(UnFixed.begin()+i);
                break;
            }
        }
    }

    s.ofv = Decoder(s, data);
}

void GRASP(const TRunData &runData, const TProblemData &data)
{
    const char* method = "GRASP";

    TSol s;
    TSol sLine;
    TSol sLineBest;
    TSol sBest;

    float currentTime = 0;
    int improv = 0;

    double start_timeMH = get_time_in_seconds();
    double end_timeMH = get_time_in_seconds();

    std::vector<int> RKorder;
    RKorder.resize(data.n);
    std::iota(RKorder.begin(), RKorder.end(), 0);

    float alphaGrasp = 0.1;
    float h;
    float hs = 0;
    float he = 0;

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

    numPar = 3;
    std::vector<std::vector<double>> parameters;
    parameters.resize(numPar);

    readParameters(method, runData.control, parameters, numPar);

    if (runData.control == 0){

        alphaGrasp = parameters[0][0];
        hs         = parameters[1][0];
        he         = parameters[2][0];
    }

    else {
        if (runData.control == 1){

            CreateStates(parameters, numStates, numPar, S);

            restartEpsilon = 1;

            epsilon_max = 1.0;

            iCurr = irandomico(0,numStates-1);

            alphaGrasp = S[iCurr].par[0];
            hs = S[iCurr].par[1];
            he = S[iCurr].par[2];
        }
    }

    CreateInitialSolutions(s, data.n);
    s.ofv = Decoder(s, data);
    sBest = s;

    int iter = 0;

    while (currentTime < runData.MAXTIME*runData.restart)
    {

        if (runData.control == 1){

            SetQLParameter(currentTime, Ti, restartEpsilon, epsilon_max, epsilon_min, epsilon, lf, df, runData.MAXTIME*runData.restart);

            at = ChooseAction(S, st, epsilon);

            iCurr = S[st].Ai[at];

            alphaGrasp = S[iCurr].par[0];
            hs = S[iCurr].par[1];
            he = S[iCurr].par[2];
        }

        h = hs;

        while (h >= he && currentTime < runData.MAXTIME*runData.restart)
        {
            if (stop_execution.load()) return;

            iter++;

            if (runData.control == 0){
                alphaGrasp = randomico(0.1, 0.9);
            }

            sLine = s;
            ConstrutiveGreedyRandomized(sLine, h, alphaGrasp, data);

            sLineBest = sLine;
            RVND(sLineBest, data, runData.strategy, RKorder);

            if (sLineBest.ofv < sBest.ofv){
                sBest = sLineBest;
                improv = 1;

                UpdatePoolSolutions(sLineBest, method, runData.debug);
            }

            else{
                h = h/2;
            }

            if (sLineBest.ofv < s.ofv){
                s = sLineBest;
            }
            else{

                if ( randomico(0.0,1.0) < (exp(-(sLineBest.ofv - s.ofv)/(100 - 100*(currentTime / runData.MAXTIME*runData.restart)))) ){
                    s = sLineBest;
                }
            }

            end_timeMH = get_time_in_seconds();
            currentTime = end_timeMH - start_timeMH;
        }

        if (runData.control == 1){

            if (improv){
                R = 1;
                improv = 0;
            }
            else{
                R = (sBest.ofv - sLineBest.ofv)/sLineBest.ofv;
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

}

#endif
