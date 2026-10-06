#ifndef _VNS_H
#define _VNS_H

void VNS(const TRunData &runData, const TProblemData &data)
{
    const char* method = "VNS";
    double beta = 0;
    int Iter = 0;
    int IterMelhora = 0;
    int kMax = 0;
    double betaMin = 0.0;
    int improv = 0;

    TSol s,
         sBest,
         sLine,
         sBestLine;

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

    numPar = 2;
    std::vector<std::vector<double>> parameters;
    parameters.resize(numPar);

    readParameters(method, runData.control, parameters, numPar);

    if (runData.control == 0){
        kMax    = parameters[0][0];
        betaMin = parameters[1][0];
    }

    else{

        if (runData.control == 1){

            CreateStates(parameters, numStates, numPar, S);

            restartEpsilon = 1;

            epsilon_max = 1.0;

            iCurr = irandomico(0,numStates-1);

            kMax    = S[iCurr].par[0];
            betaMin = S[iCurr].par[1];
        }
    }

    sBest.ofv = INFINITY;
    for (int i=0; i<1; i++)
    {
        CreateInitialSolutions(s, data.n);
        s.ofv = Decoder(s, data);
        if (s.ofv < sBest.ofv)
            sBest = s;
    }

    s = sBest;

    while (currentTime < runData.MAXTIME*runData.restart)
    {

        if (runData.control == 1){

            SetQLParameter(currentTime, Ti, restartEpsilon, epsilon_max, epsilon_min, epsilon, lf, df, runData.MAXTIME*runData.restart);

            at = ChooseAction(S, st, epsilon);

            iCurr = S[st].Ai[at];

            kMax    = S[iCurr].par[0];
            betaMin = S[iCurr].par[1];
        }

        int k = 1;
        while (k <= kMax && currentTime < runData.MAXTIME*runData.restart)
        {
            if (stop_execution.load()) return;

            Iter++;

            beta = randomico(k*betaMin,(k+1)*betaMin);

            sLine = s;
            ShakeSolution(sLine, beta, beta, data.n);

            sLine.ofv = Decoder(sLine, data);

            sBestLine = sLine;
            RVND(sBestLine, data, runData.strategy, RKorder);

            if (sBestLine.ofv < s.ofv)
            {
                s = sBestLine;

                if (s.ofv < sBest.ofv){
                    sBest = s;

                    k = 1;
                    IterMelhora = Iter;
                    improv = 1;

                    UpdatePoolSolutions(sBestLine, method, runData.debug);
                }
            }
            else
            {

                k++;

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
                R = (sBest.ofv - s.ofv)/s.ofv;
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
    }

}

#endif
