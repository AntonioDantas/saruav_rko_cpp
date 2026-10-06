#ifndef _ILS_H
#define _ILS_H

void ILS(const TRunData &runData, const TProblemData &data)
{
    const char* method = "ILS";
    int Iter = 0;
    int IterImprov = 0;

    TSol sBest,
         sLine,
         sBestLine;

    double betaMin = 0.0;
    double betaMax = 0.0;

    float currentTime = 0;
    int improv = 0;

    double start_timeMH = get_time_in_seconds();
    double end_timeMH = get_time_in_seconds();

    std::vector<int> RKorder;
    RKorder.resize(data.n);
    for (int i = 0; i < data.n; i++){
        RKorder[i] = i;
    }

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
        betaMin = parameters[0][0];
        betaMax = parameters[1][0];
    }

    else{

        if (runData.control == 1){

            CreateStates(parameters, numStates, numPar, S);

            restartEpsilon = 1;

            epsilon_max = 1.0;

            iCurr = irandomico(0,numStates-1);

            betaMin = S[iCurr].par[0];
            betaMax = S[iCurr].par[1];
        }
    }

    Iter = 0;

    CreateInitialSolutions(sBest, data.n);
    sBest.ofv = Decoder(sBest, data);

    RVND(sBest, data, runData.strategy, RKorder);
    UpdatePoolSolutions(sBest, method, runData.debug);

    end_timeMH = get_time_in_seconds();
    currentTime = end_timeMH - start_timeMH;

    while (currentTime < runData.MAXTIME*runData.restart)
    {
        if (stop_execution.load()) return;

        Iter++;

        if (runData.control == 1){

            SetQLParameter(currentTime, Ti, restartEpsilon, epsilon_max, epsilon_min, epsilon, lf, df, runData.MAXTIME*runData.restart);

            at = ChooseAction(S, st, epsilon);

            iCurr = S[st].Ai[at];

            betaMin = S[iCurr].par[0];
            betaMax = S[iCurr].par[1];
        }

        sLine = sBest;

        ShakeSolution(sLine, betaMin, betaMax, data.n);

        sLine.ofv = Decoder(sLine, data);

        sBestLine = sLine;
        RVND(sBestLine, data, runData.strategy, RKorder);

        if (sBestLine.ofv < sBest.ofv)
        {
            sBest = sBestLine;
            IterImprov = Iter;
            improv = 1;

            UpdatePoolSolutions(sBest, method, runData.debug);
        }

        if (runData.control == 1){

            if (improv){
                R = 1;
                improv = 0;
            }
            else{
                R = (sBest.ofv - sBestLine.ofv)/sBestLine.ofv;
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
