#ifndef _QLearning_H
#define _QLearning_H

void PrintPolicy(std::vector<TState> S, int st)
{
    printf("\n\nPolicy\n");
    for (int i=0; i<(int)S.size(); i++)
    {
        printf("\n i = %d: \t", i);
        for (int j=0; j<(int)S[i].Ai.size(); j++)
        {
            printf("%d (%.2lf) \t", S[i].Ai[j], S[i].Qa[j]);
        }
    }
}

void CreateStates(std::vector<std::vector<double>> parameters, int &numStates, int numPar, std::vector<TState> &S)
{
    #pragma omp critical
    {

        numStates = parameters[0].size();
        for (int i=1; i<numPar; i++){
            numStates = numStates * parameters[i].size();
        }

        for (int i=0; i<numStates; i++)
        {
            TState sAux;
            sAux.label = i;
            sAux.par.resize(numPar);
            sAux.ci = 0;
            sAux.numN = 0;
            sAux.Ai.clear();
            sAux.Qa.clear();
            sAux.sumQ = 0;
            sAux.maxQ = 0;
            sAux.minA = 0;
            sAux.maxA = 0;

            S.push_back(sAux);
        }

        int nC = 1;
        for (int j=0; j<numPar; j++)
        {
            nC = nC * parameters[j].size();

            for (int i=0; i<numStates; i++)
            {
                int columnIndex = i / (numStates/nC);
                columnIndex = columnIndex % parameters[j].size();

                S[i].par[j] = parameters[j][columnIndex];
            }
        }

        for (int i=0; i<(int)S.size(); i++)
        {
            for (int j=0; j<(int)S.size(); j++)
            {

                int distance = 0;
                for (int k=0; k<(int)S[i].par.size(); k++) {
                    if (S[i].par[k] != S[j].par[k]) {
                        distance++;
                    }
                }

                if (distance <= 1){
                    S[i].Ai.push_back(j);

                    int q0 = randomico(0.05,0.01);
                    S[i].Qa.push_back(q0);

                    if (q0 > S[i].maxQ)
                    {
                        S[i].maxQ = q0;
                        S[i].maxA = S[i].Ai.size()-1;
                    }
                }
            }
        }

    }
}

int ChooseAction(std::vector<TState> S, int st, double epsilon)
{

    int at;

    if (randomico(0,1) <= 1-epsilon)
    {

        at = S[st].maxA;
    }
    else
    {

        at = irandomico(0,S[st].Ai.size()-1);
    }

    return at;
}

void SetQLParameter(float currentTime, int &Ti, int &restartEpsilon, float epsilon_max, float epsilon_min, double &epsilon, double &lf, double &df, int MAXTIME)
{

    static const double PI = 3.14159265;

    Ti = MAXTIME * 0.1;
    if (currentTime >= restartEpsilon * Ti){
        restartEpsilon++;

        epsilon_max = epsilon_max - 0.1;
        if (epsilon_max < epsilon_min)
            epsilon_max = epsilon_min;
        epsilon = epsilon_max;
    }
    else {
        epsilon = epsilon_min + 0.5 * (epsilon_max - epsilon_min) * (1 + cos((((int)currentTime%Ti)/(float)(Ti))*PI));
    }

    lf = 1 - (0.9 * currentTime / MAXTIME);

    df = 0.8;
}

#endif
