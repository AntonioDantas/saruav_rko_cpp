#ifndef _GA_H
#define _GA_H

static void UpdatePopulation(int sizePop, std::vector <TSol> &Pop, const TProblemData &data)
{

    int oldPsize = Pop.size();

    if (oldPsize > sizePop){

        Pop.resize(sizePop);
    }

    else if (oldPsize < sizePop){

        Pop.resize(sizePop);

        for (int k = oldPsize; k < sizePop; k++)
        {
            CreateInitialSolutions(Pop[k], data.n);
            Pop[k].ofv = Decoder(Pop[k], data);
        }
    }
}

void GA(const TRunData &runData, const TProblemData &data)
{
    const char* method = "GA";
    int sizePop = 0;
    float probCros = 0.0;
    float probMut = 0.0;

    std::vector <TSol> Pop;
    std::vector <TSol> PopInter;
    std::vector <TSol> PopNew;

    TSol bestInd;

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

    numPar = 3;
    std::vector<std::vector<double>> parameters;
    parameters.resize(numPar);

    readParameters(method, runData.control, parameters, numPar);

    if (runData.control == 0){

        sizePop  = parameters[0][0];
        probCros = parameters[1][0];
        probMut  = parameters[2][0];
    }

    else {

        if (runData.control == 1){

            CreateStates(parameters, numStates, numPar, S);

            restartEpsilon = 1;

            epsilon_max = 1.0;

            iCurr = irandomico(0,numStates-1);

            sizePop  = (int)S[iCurr].par[0];
            probCros = S[iCurr].par[1];
            probMut  = S[iCurr].par[2];
        }
    }

    Pop.clear();
    PopNew.clear();
    PopInter.clear();

    Pop.resize(sizePop);
    PopNew.resize(sizePop);
    PopInter.resize(sizePop);

    bestInd.ofv = 9999999999;
    for (int i=0; i<sizePop; i++)
    {
        CreateInitialSolutions(Pop[i], data.n);
        Pop[i].ofv = Decoder(Pop[i], data);

        if (Pop[i].ofv < bestInd.ofv)
            bestInd = Pop[i];
    }

    while (currentTime < runData.MAXTIME*runData.restart)
    {

        numGenerations++;

        if (runData.control == 1){

            SetQLParameter(currentTime, Ti, restartEpsilon, epsilon_max, epsilon_min, epsilon, lf, df, runData.MAXTIME*runData.restart);

            at = ChooseAction(S, st, epsilon);

            iCurr = S[st].Ai[at];

            sizePop  = (int)S[iCurr].par[0];
            probCros = S[iCurr].par[1];
            probMut  = S[iCurr].par[2];

            UpdatePopulation(sizePop, Pop, data);
        }

        PopInter.clear();
        PopInter.resize(sizePop);
        for (int j=0; j<sizePop; j++)
        {
            int p1, p2, p3, melhor;

            do{
                p1 = irandomico(0,sizePop-1);
                p2 = irandomico(0,sizePop-1);
                p3 = irandomico(0,sizePop-1);
            }while (p1 == p2 || p1 == p3 || p2 == p3);

            if ((Pop[p1].ofv < Pop[p2].ofv) && (Pop[p1].ofv < Pop[p3].ofv))
                melhor = p1;
            else
            if (Pop[p2].ofv < Pop[p3].ofv)
                melhor = p2;
            else
                melhor = p3;

            PopInter[j] = Pop[melhor];
        }

        double avrgOff = 0.0;
        double bestOFcurrent = INFINITY;
        PopNew.clear();
        PopNew.resize(sizePop);
        for (int i=0; i<sizePop-1; i=i+2)
        {
            if (stop_execution.load()) return;

            PopNew[i] = PopInter[i];
            PopNew[i+1] = PopInter[i+1];

            if (randomico(0,1) < probCros)
            {
                for(int j = 0; j < data.n; j++)
                {

                    if (randomico(0,1) < 0.5){
                        PopNew[i].rk[j] = PopInter[i+1].rk[j];
                    }

                    if (randomico(0,1) < 0.5){
                        PopNew[i+1].rk[j] = PopInter[i].rk[j];
                    }

                    if (randomico(0,1) <= probMut){
                        PopNew[i].rk[j] = randomico(0,1);
                    }

                    if (randomico(0,1) <= probMut){
                        PopNew[i+1].rk[j] = randomico(0,1);
                    }
                }

                PopNew[i].ofv   = Decoder(PopNew[i], data);
                PopNew[i+1].ofv = Decoder(PopNew[i+1], data);
            }

            avrgOff += PopNew[i].ofv + PopNew[i+1].ofv;

            if (PopNew[i].ofv < bestInd.ofv){
                bestInd = PopNew[i];
                bestGeneration = numGenerations;
                improv = 1;

                UpdatePoolSolutions(bestInd, method, runData.debug);
            }

            if (PopNew[i+1].ofv < bestInd.ofv){
                bestInd = PopNew[i+1];
                bestGeneration = numGenerations;
                improv = 1;

                UpdatePoolSolutions(bestInd, method, runData.debug);
            }

            if (PopNew[i].ofv < bestOFcurrent){
                bestOFcurrent = PopNew[i].ofv;
            }
            if (PopNew[i+1].ofv < bestOFcurrent){
                bestOFcurrent = PopNew[i+1].ofv;
            }
        }

        int pos1 = irandomico(0,sizePop-1);
        NelderMeadSearch(PopNew[pos1], data);

        if (PopNew[pos1].ofv < bestInd.ofv){
            bestInd = PopNew[pos1];
            bestGeneration = numGenerations;

            UpdatePoolSolutions(bestInd, method, runData.debug);
        }

        Pop = PopNew;

        if (runData.control == 1){

            if (improv){
                R = 1;
                improv = 0;
            }
            else{
                R = (bestInd.ofv - bestOFcurrent)/bestOFcurrent;
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
    PopNew.clear();
    PopInter.clear();

}

#endif
