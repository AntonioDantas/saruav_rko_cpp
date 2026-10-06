#ifndef _BRKGA_CS_H
#define _BRKGA_CS_H

void BRKGA_CS(const TRunData &runData, const TProblemData &data);

static void UpdatePopulationSize(int p, double pe, double pm, double rhoe, std::vector <TSol> &Pop, std::vector <TSol> &PopInter, const TProblemData &data);

static void ChaoticInd(TSol &s, int rhoe, const int n);

static TSol ParametricUniformCrossover(int elitesize, int popSize, double pm, double rhoe, std::vector <TSol> Pop, const int n);

static double PearsonCorrelation(std::vector <double> s1, std::vector <double> s2, const int n);

static void IC(int p, double pe, std::vector <TSol> &Pop, std::vector<int> &promising, const int n);

static void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas, std::vector <TSol> &Pop, std::vector<int> &label);

static void PromisingLP(int p, double pe, std::vector <TSol> &Pop, std::vector<int> &label, std::vector<int> &promising);

void BRKGA_CS(const TRunData &runData, const TProblemData &data)
{
    const char* method = "BRKGA-CS";

    int p = 1597;
    double pe = 0.20;
    double pm = 0.03;
    double rhoe = 0.70;

    std::vector <TSol> Pop;
    std::vector <TSol> PopInter;
    TSol bestInd;

    double start_timeMH = get_time_in_seconds();
    double end_timeMH = get_time_in_seconds();

    int numGenerations = 0;
    int bestGeneration = 0;
    float currentTime = 0;
    int sumLS = 0;
    int noImprovBRKGA = 0;
    int improv = 0;

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

    if (runData.control == 0)
    {
        p    = parameters[0][0];
        pe   = parameters[1][0];
        pm   = parameters[2][0];
        rhoe = parameters[3][0];
    }

    else
    {

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

        noImprovBRKGA++;

        if (runData.control == 1){

            SetQLParameter(currentTime, Ti, restartEpsilon, epsilon_max, epsilon_min, epsilon, lf, df, runData.MAXTIME*runData.restart);

            at = ChooseAction(S, st, epsilon);

            iCurr = S[st].Ai[at];

            p       = (int)S[iCurr].par[0];
            pe      = S[iCurr].par[1];
            pm      = S[iCurr].par[2];
            rhoe    = S[iCurr].par[3];

            UpdatePopulationSize(p, pe, pm, rhoe, Pop, PopInter, data);
        }

        for (int i=0; i<(int)(p*pe); i++){

            PopInter[i] = Pop[i];
        }

        double bestOFV = INFINITY;
        for (int i = (int)(p*pe); i < p; i++){
            if (stop_execution.load()) return;

            PopInter[i] = ParametricUniformCrossover((int)(p*pe), p, pm, rhoe, Pop, data.n);

            PopInter[i].ofv = Decoder(PopInter[i], data);

            if (PopInter[i].ofv < bestOFV)
                bestOFV = PopInter[i].ofv;
        }

        Pop = PopInter;

        sort(Pop.begin(), Pop.end(), sortByFitness);

        if (Pop[0].ofv < bestInd.ofv){
            bestInd = Pop[0];
            bestGeneration = numGenerations;
            noImprovBRKGA = 0;

            UpdatePoolSolutions(bestInd, method, runData.debug);
        }

        if (runData.control == 1){

            if (improv){

                R = 1 + 1/p;
                improv = 0;
            }
            else{
                R = (bestInd.ofv - bestOFV)/bestOFV;
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

        sumLS = 0;

        if (noImprovBRKGA == 0 || noImprovBRKGA > 10)
        {

            noImprovBRKGA = 0;

            std::vector<int> promising((int)(p*pe), 0);

            IC(p, pe, Pop, promising, data.n);

            std::vector <int> promisingSol;
            promisingSol.clear();

            for (int i=0; i < (int)(p*pe); i++) {
                if (stop_execution.load()) return;

                if (promising[i] == 1){
                    promisingSol.push_back(i);
                }

                else {
                    ChaoticInd(Pop[i], rhoe, data.n);
                    Pop[i].ofv = Decoder(Pop[i], data);
                }
            }

            for (unsigned int i=0; i < promisingSol.size(); i++){
                if (stop_execution.load()) return;

                if (i < 1)
                   RVND(Pop[promisingSol[i]], data, runData.strategy, RKorder);
                else
                   NelderMeadSearch(Pop[promisingSol[i]], data);

                if (Pop[promisingSol[i]].ofv < bestInd.ofv){
                    bestInd = Pop[promisingSol[i]];
                    bestGeneration = numGenerations;
                    noImprovBRKGA = 0;

                    UpdatePoolSolutions(bestInd, method, runData.debug);
                }

            }

            sumLS = promisingSol.size();
            promisingSol.clear();
            promising.clear();

            sort(Pop.begin(), Pop.end(), sortByFitness);
        }

        end_timeMH = get_time_in_seconds();
        currentTime = end_timeMH - start_timeMH;
    }

    Pop.clear();
    PopInter.clear();
}

static void UpdatePopulationSize(int p, double pe, double pm, double rhoe, std::vector <TSol> &Pop, std::vector <TSol> &PopInter, const TProblemData &data)
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

        	Pop[k] = ParametricUniformCrossover((int)(oldPsize*pe), oldPsize-1, pm, rhoe, Pop, data.n);
            Pop[k].ofv = Decoder(Pop[k], data);
        }

        sort(Pop.begin(), Pop.end(), sortByFitness);

        PopInter.clear();
        PopInter.resize(p);
    }
}

static void ChaoticInd(TSol &s, int rhoe, const int n)
{

    for (int k=0; k<n; k++)
    {
        if (randomico(0,1) > rhoe)
           s.rk[k] = randomico(0,1);
    }
}

static TSol ParametricUniformCrossover(int eliteSize, int popSize, double pm, double rhoe, std::vector <TSol> Pop, const int n)
{
	TSol s;

    int eliteParent = irandomico(0, eliteSize - 1);

    int nonEliteParent = irandomico(eliteSize, popSize-1);

	s = Pop[eliteParent];

    for(int j = 0; j < n; j++)
    {

        if (randomico(0,1) < pm)
        {
            s.rk[j] = randomico(0,1);
        }

        else
        {

            if (randomico(0,1) < rhoe){
                s.rk[j] = Pop[eliteParent].rk[j];
            }
            else{
                s.rk[j] = Pop[nonEliteParent].rk[j];
            }
        }
    }

    return s;
}

static double PearsonCorrelation(std::vector <double> X, std::vector <double> Y, const int n)
{
    double correlation = 0;
    double sumXY = 0;
    double sumX2 = 0;
    double sumY2 = 0;
    double sumX = 0;
    double sumY = 0;

    for(int j=0; j<n; j++)
    {
        sumX += X[j];
        sumX2 += X[j] * X[j];
        sumXY += X[j] * Y[j];
        sumY += Y[j];
        sumY2 += Y[j] * Y[j];
    }

    correlation= ((n*sumXY) - (sumX*sumY) ) / (sqrt( (n*sumX2 - sumX*sumX) * (n*sumY2 - sumY*sumY) ));
    return correlation;
}

static void IC(int p, double pe, std::vector <TSol> &Pop, std::vector<int> &promising, const int n)
{
    int Tpe = (int)p*pe;
    std::vector<std::vector<std::pair<int, double> > > listaArestas(Tpe, std::vector<std::pair<int, double> >());

    double sigma = 0.6;

	int entrouAresta = 0;
	double pearson = 0.0;
	for (int i = 0; i < Tpe - 1; i++) {
		for (int j = i + 1; j < Tpe; j++)
		{
			pearson = PearsonCorrelation(Pop[i].rk, Pop[j].rk, n);
			if (pearson > sigma) {
				entrouAresta++;
				listaArestas[i].push_back(std::make_pair(j, pearson));
				listaArestas[j].push_back(std::make_pair(i, pearson));
			}
            else{
                entrouAresta += 5;
            }
		}
	}

    std::vector<int> label(Tpe,0);

	LP(listaArestas, Pop, label);

	PromisingLP(p, pe, Pop, label, promising);
    listaArestas.clear();
    label.clear();
}

static void LP(std::vector<std::vector<std::pair<int, double> > > listaArestas, std::vector <TSol> &Pop, std::vector<int> &label)
{
    int nk = listaArestas.size();

	std::vector<int> ordemVisita(nk);
	iota(ordemVisita.begin(), ordemVisita.end(), 0);

	for (int i = 0; i < nk; i++)
		label[i] = i;

	int iteracao = 1;
	int labelVizinho, melhorLabel;
	double melhorPeso;
	std::map<int, double> totalLabels;
	std::map<int, double>::iterator it;

	int movimentos = 1;
	while (movimentos)
    {
		movimentos = 0;
		shuffle(ordemVisita.begin(), ordemVisita.end(),std::mt19937(std::random_device()()));
        for (std::vector<int>::size_type idVertice=0; idVertice <ordemVisita.size(); idVertice++)
        {

			totalLabels.clear();
            for (std::vector<std::pair<int, double> >::iterator itVizinho = listaArestas[idVertice].begin();
                itVizinho != listaArestas[idVertice].end(); ++itVizinho) {
                int idVizinho = itVizinho->first;
                labelVizinho = label[idVizinho];
                it = totalLabels.find(labelVizinho);
                if (it != totalLabels.end()) {
                    it->second += itVizinho->second;
                } else {
                    totalLabels[labelVizinho] = itVizinho->second;
                }
            }

            melhorLabel = label[idVertice];
			melhorPeso = std::numeric_limits<double>::min();
            for (std::map<int, double>::iterator itTotais = totalLabels.begin(); itTotais != totalLabels.end(); ++itTotais) {
                if (itTotais->second > melhorPeso) {
                    melhorLabel = itTotais->first;
                    melhorPeso = itTotais->second;
                }
            }

			if (melhorLabel != label[idVertice]) {
				label[idVertice] = melhorLabel;
				movimentos = 1;
			}
		}
		iteracao++;
	}

    ordemVisita.clear();
}

static void PromisingLP(int p, double pe, std::vector <TSol> &Pop, std::vector<int> &label, std::vector<int> &promising)
{
    int Tpe = (int)p*pe;
    std::vector <int> grupos;
	int tamanhoGrupos = 0;

	int found;

    for (int i = 0; i < Tpe; i++)
	{
		found = 0;
		for (unsigned int j = 0; j < grupos.size(); j++)
		{
			if (label[i] == grupos[j])
                found = 1;
		}
		if (found == 0)
		{
			tamanhoGrupos++;
			grupos.push_back(label[i]);
		}
	}

    for (unsigned int j = 0; j < grupos.size(); j++)
	{
        int numTry = 0;
        int pos = 0;

        while(numTry < Tpe)
        {
            numTry++;
            pos = irandomico(0,Tpe-1);

            if (label[pos] == grupos[j])
            {
                promising[pos] = 1;
                break;
			}
        }
    }
}

#endif
