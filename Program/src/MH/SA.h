#ifndef _SA_H
#define _SA_H

void SA(const TRunData &runData, const TProblemData &data)
{
    const char* method = "SA";
    TSol s;
    TSol sViz;
    TSol sBest;

    double delta = 0;
    double bestOFV;

    double T0 = 0;
    double T = 0;
    double alphaSA = 0;
    int SAmax = 0;
    int IterT = 0;
    float betaMin = 0;
    float betaMax = 0;
    int reanneling = 0;
    int improv = 0;
    float currentTime = 0;

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

    numPar = 5;
    std::vector<std::vector<double>> parameters;
    parameters.resize(numPar);

    readParameters(method, runData.control, parameters, numPar);

    if (runData.control == 0){

        SAmax   = parameters[0][0];
        alphaSA = parameters[1][0];
        betaMin = parameters[2][0];
        betaMax = parameters[3][0];
        T0      = parameters[4][0];
    }

    else
    {

        if (runData.control == 1){

            CreateStates(parameters, numStates, numPar, S);

            restartEpsilon = 1;

            epsilon_max = 1.0;

            iCurr = irandomico(0,numStates-1);

            T0 = 1000000;
            SAmax = (int)S[iCurr].par[0];
            alphaSA = S[iCurr].par[1];
            betaMin = S[iCurr].par[2];
            betaMax = S[iCurr].par[3];
            T0      = S[iCurr].par[4];
        }
    }

    CreateInitialSolutions(s, data.n);
    s.ofv = Decoder(s, data);
    sBest = s;

    while (currentTime < runData.MAXTIME*runData.restart)
    {
        IterT = 0;
        if (!reanneling) T = T0;
        else T = T0*0.3;
        while (T > 0.0001 && currentTime < runData.MAXTIME*runData.restart)
        {

            if (runData.control == 1){

                SetQLParameter(currentTime, Ti, restartEpsilon, epsilon_max, epsilon_min, epsilon, lf, df, runData.MAXTIME*runData.restart);

                at = ChooseAction(S, st, epsilon);

                iCurr = S[st].Ai[at];

                SAmax   = (int)S[iCurr].par[0];
                alphaSA = S[iCurr].par[1];
                betaMin = S[iCurr].par[2];
                betaMax = S[iCurr].par[3];
            }

            bestOFV = INFINITY;
            while (IterT < SAmax && currentTime < runData.MAXTIME*runData.restart)
            {
                if (stop_execution.load()) return;

                IterT++;

                sViz = s;
                ShakeSolution(sViz, betaMin, betaMax, data.n);

                sViz.ofv = Decoder(sViz, data);

                if (sViz.ofv < bestOFV)
                    bestOFV = sViz.ofv;

                delta = sViz.ofv - s.ofv;

                if (delta < 0)
                {

                    s = sViz;

                    if (s.ofv < sBest.ofv)
                    {
                        sBest = s;
                        improv = 1;

                        UpdatePoolSolutions(s, method, runData.debug);
                    }
                }
                else
                {

                    double x = randomico(0,1);

                    if ( x < (exp(-delta/T)) )
                        s = sViz;
                }
            }

            if (runData.control == 1){

                if (improv){
                    R = 1;
                    improv = 0;
                }
                else{
                    R = (sBest.ofv - bestOFV)/bestOFV;
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

            T = T * alphaSA;
            IterT = 0;

            sViz = s;

            NelderMeadSearch(sViz, data);

            if (sViz.ofv < sBest.ofv)
            {
                sBest = sViz;

                UpdatePoolSolutions(sBest, method, runData.debug);
            }

            end_timeMH = get_time_in_seconds();
            currentTime = end_timeMH - start_timeMH;

        }

        reanneling = 1;
    }

}

#endif
