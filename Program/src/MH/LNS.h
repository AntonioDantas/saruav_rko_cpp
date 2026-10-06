#ifndef _LNS_H
#define _LNS_H

void fareySequence(int num, std::vector<double> &F) {
    int a = 0, b = 1, c = 1, d = num;

    F.push_back((double)a/b);

    while (c <= num) {
        int k = (num + b) / d;
        int temp_a = a, temp_b = b;
        a = c;
        b = d;
        c = k * c - temp_a;
        d = k * d - temp_b;
        if (a != b){

            F.push_back((double)a/b);
        }
        else{
            F.push_back(((double)a/b) - 0.00001);
        }
    }
}

void LNS(const TRunData &runData, const TProblemData &data)
{
    const char* method = "LNS";
    double T0 = 0;
    double T;
    double alphaLNS = 0;
    double betaMin = 0.0;
    double betaMax = 0.0;
    int reanneling = 0;

    TSol sLine,
         sLineBest,
         s,
         sBest;

    float currentTime = 0;
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

    std::vector<double> F;
    fareySequence(7,F);

    if (runData.control == 0){
        betaMin  = parameters[0][0];
        betaMax  = parameters[1][0];
        T0       = parameters[2][0];
        alphaLNS = parameters[3][0];
    }

    else{

        if (runData.control == 1){

            CreateStates(parameters, numStates, numPar, S);

            restartEpsilon = 1;

            epsilon_max = 1.0;

            iCurr = irandomico(0,numStates-1);

            betaMin  = S[iCurr].par[0];
            betaMax  = S[iCurr].par[1];
            T0       = S[iCurr].par[2];
            alphaLNS = S[iCurr].par[3];
        }
    }

    CreateInitialSolutions(s, data.n);
    s.ofv = Decoder(s, data);
    sBest = s;

    while (currentTime < runData.MAXTIME*runData.restart)
    {
        if (!reanneling) T = T0;
        else T = T0*0.3;
        while (T > 0.01 && currentTime < runData.MAXTIME*runData.restart)
        {

            if (runData.control == 1){

                SetQLParameter(currentTime, Ti, restartEpsilon, epsilon_max, epsilon_min, epsilon, lf, df, runData.MAXTIME*runData.restart);

                at = ChooseAction(S, st, epsilon);

                iCurr = S[st].Ai[at];

                betaMin  = S[iCurr].par[0];
                betaMax  = S[iCurr].par[1];
                T0       = S[iCurr].par[2];
                alphaLNS = S[iCurr].par[3];
            }

            sLine = s;
            int intensity = irandomico(betaMin*data.n, betaMax*data.n);

            std::shuffle (RKorder.begin(), RKorder.end(),rng);

            for (int k=0; k<intensity; k++)
            {
                int pos = RKorder[k];
                double OFVbest = INFINITY;
                double rkBest = 0;

                for (int j=0; j<(int)F.size()-1; j++)
                {
                    if (stop_execution.load()) return;

                    sLine.rk[pos] = randomico(F[j], F[j+1]);
                    sLine.ofv = Decoder(sLine, data);

                    if (sLine.ofv < OFVbest)
                    {
                        OFVbest = sLine.ofv;
                        rkBest = sLine.rk[pos];
                    }
                }

                sLine.ofv = OFVbest;
                sLine.rk[pos] = rkBest;
            }

            sLineBest = sLine;
            NelderMeadSearch(sLineBest, data);

            double delta = sLineBest.ofv - s.ofv;

            if (delta <= 0)
            {
                s = sLineBest;

                if (s.ofv < sBest.ofv){
                    sBest = s;
                    improv = 1;

                    UpdatePoolSolutions(s, method, runData.debug);
                }
            }
            else
            {
                double x = randomico(0,1);

                if ( x < (exp(-delta/T)) )
                    s = sLineBest;
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

            T = T * alphaLNS;

            end_timeMH = get_time_in_seconds();
            currentTime = end_timeMH - start_timeMH;

        }

        reanneling = 1;
    }

}

#endif
