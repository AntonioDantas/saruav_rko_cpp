#ifndef _BRKGA_H
#define _BRKGA_H

static TSol PUX(int eliteSize, int popSize, double rhoe, std::vector <TSol> &Pop, const int n)
{
	TSol s;

    int eliteParent = irandomico(0, eliteSize - 1);
    int nonEliteParent = irandomico(eliteSize, popSize-1);

	s.rk.resize(n);

    for(int j = 0; j < n; j++)
    {

        if (randomico(0,1) < rhoe){
            s.rk[j] = Pop[eliteParent].rk[j];
        }
        else{
            s.rk[j] = Pop[nonEliteParent].rk[j];
        }
    }

    return s;
}

static void UpdatePopSize(int p, double pe, double pm, double rhoe, std::vector <TSol> &Pop, std::vector <TSol> &PopInter, const TProblemData &data)
{

    int oldPsize = Pop.size();

    if (oldPsize > p){

        PopInter = Pop;

        Pop.resize(p);

        for (int i=0; i<(int)(p*pe); i++){

            Pop[i] = PopInter[i];
        }

        int pos = (int)(pe*oldPsize);
        for (int i=(int)(p*pe); i<p; i++){

            Pop[i] = PopInter[pos];
            pos++;
        }

        PopInter.clear();
        PopInter.resize(p);
    }

    else if (oldPsize < p){

        Pop.resize(p);

        for (int k = oldPsize; k < p; k++)
        {
            if (stop_execution.load()) return;

        	Pop[k] = PUX((int)(oldPsize*pe), oldPsize-1, rhoe, Pop, data.n);
            Pop[k].ofv = Decoder(Pop[k], data);
        }

        sort(Pop.begin(), Pop.end(), sortByFitness);

        PopInter.clear();
        PopInter.resize(p);
    }
}

void BRKGA(const TRunData &runData, const TProblemData &data)
{
    const char* method = "BRKGA";
    int p = 0;
    double pe = 0.0;
    double pm = 0.0;
    double rhoe = 0.0;

    std::vector <TSol> Pop;
    std::vector <TSol> PopInter;
    TSol bestInd;

    int numGenerations = 0;
    int bestGeneration = 0;
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

    if (runData.control == 0){
        p    = parameters[0][0];
        pe   = parameters[1][0];
        pm   = parameters[2][0];
        rhoe = parameters[3][0];
    }

    else {

        if (runData.control == 1){

            CreateStates(parameters, numStates, numPar, S);

            restartEpsilon = 1;

            epsilon_max = 1.0;

            iCurr = irandomico(0,numStates-1);

            p    = (int)S[iCurr].par[0];
            pe   = S[iCurr].par[1];
            pm   = S[iCurr].par[2];
            rhoe = S[iCurr].par[3];
        }
    }

    Pop.clear();
    PopInter.clear();
    Pop.resize(p);
    PopInter.resize(p);

    for (int i=0; i<p; i++)
    {
        CreateInitialSolutions(Pop[i], data.n);
        Pop[i].ofv = Decoder(Pop[i], data);
        PopInter[i] = Pop[i];
    }

    sort(Pop.begin(), Pop.end(), sortByFitness);

    bestInd = Pop[0];

    while (currentTime < runData.MAXTIME*runData.restart)
    {

        numGenerations++;

        if (runData.control == 1){

            SetQLParameter(currentTime, Ti, restartEpsilon, epsilon_max, epsilon_min, epsilon, lf, df, runData.MAXTIME*runData.restart);

            at = ChooseAction(S, st, epsilon);

            iCurr = S[st].Ai[at];

            p       = (int)S[iCurr].par[0];
            pe      = S[iCurr].par[1];
            pm      = S[iCurr].par[2];
            rhoe    = S[iCurr].par[3];

            UpdatePopSize(p, pe, pm, rhoe, Pop, PopInter, data);
        }

        for (int i=0; i<(int)(p*pe); i++){

            PopInter[i] = Pop[i];
        }

        TSol bestOff;
        bestOff.ofv = INFINITY;
        for (int i = (int)(p*pe); i < p - (int)(p*pm); i++){
            if (stop_execution.load()) return;

            PopInter[i] = PUX((int)(p*pe), p, rhoe, Pop, data.n);

            PopInter[i].ofv = Decoder(PopInter[i], data);

            if (PopInter[i].ofv < bestOff.ofv){
                bestOff = PopInter[i];
            }
        }

        for (int i = p - (int)(p*pm) - (int)(p*pe); i < p; i++){
            if (stop_execution.load()) return;

            CreateInitialSolutions(PopInter[i], data.n);
            PopInter[i].ofv = Decoder(PopInter[i], data);

        }

        Pop = PopInter;

        int pos = irandomico(0, (int)p*pe);
        NelderMeadSearch(Pop[pos], data);

        sort(Pop.begin(), Pop.end(), sortByFitness);

        if (Pop[0].ofv < bestInd.ofv){
            bestGeneration = numGenerations;
            improv = 1;
            bestInd = Pop[0];

            UpdatePoolSolutions(Pop[0], method, runData.debug);
        }

        if (runData.control == 1){

            if (improv){

                R = 1 + 1/p;
                improv = 0;
            }
            else{
                R = (bestInd.ofv - bestOff.ofv)/bestInd.ofv;
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

    Pop.clear();
    PopInter.clear();

}

#endif
