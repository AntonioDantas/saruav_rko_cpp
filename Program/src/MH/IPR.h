#ifndef _IPR_H
#define _IPR_H

void IPR(const TRunData &runData, const TProblemData &data)
{
    const char* method = "IPR";

    float currentTime = 0;
    double start_timeMH = get_time_in_seconds();
    double end_timeMH = get_time_in_seconds();

    TSol atual;
    TSol guia;
    TSol bestPath;
    TSol bestIteration;
    TSol sCurrent;
    TSol sViz;

    bestPath.ofv = std::numeric_limits<double>::infinity();
    bestPath.best_time = 0.0;

    while (currentTime < runData.MAXTIME*runData.restart)
    {

        int k1, k2;

        #pragma omp critical
        {
            do {
                k1 = irandomico(0,pool.size()-1);
                k2 = irandomico(0,pool.size()-1);
            }
            while (k1 == k2);

            atual = pool[k1];
            guia = pool[k2];
        }

        bestPath = atual;
        bestIteration = atual;
        sCurrent = atual;
        sViz = atual;

        int direction = 1;

        int blockSize = data.n*0.10;
        if (blockSize < 1) blockSize = 1;

        int numBlock = data.n/blockSize;
        std:: vector<int> fixedBlock(numBlock,1);
        int dist=0;

        for (int i=0; i<data.n; i++)
        {
            if (atual.rk[i] != guia.rk[i]){
                fixedBlock[i%numBlock] = 0;
                dist++;
            }
        }

        if (dist > 0)
        {

            bestPath.ofv = INFINITY;
            int bestBlock=-1;

            int numIteration = 0;
            while(numIteration < numBlock-1)
            {
                numIteration++;
                bestBlock=-1;
                bestIteration.ofv = INFINITY;

                for(int i=0; i<numBlock; i++)
                {
                    if (fixedBlock[i] == 0){
                        int initialBlock = i*blockSize;
                        int finalBlock = initialBlock + blockSize;

                        sViz = sCurrent;

                        for (int k=initialBlock; k<finalBlock && k<data.n; k++)
                        {

                            if (direction == 1){
                                sViz.rk[k] = guia.rk[k];
                            }

                            else
                            if (direction == -1){
                                sViz.rk[k] = std::max(0.0, std::min(1.0 - guia.rk[k], 0.999999));
                            }
                        }

                        sViz.ofv = Decoder(sViz, data);

                        if(sViz.ofv < bestIteration.ofv){
                            bestIteration = sViz;
                            bestBlock = i;
                        }
                    }
                }

                if (bestIteration.ofv < bestPath.ofv)
                {
                    bestPath = bestIteration;
                }

                if (bestBlock >= 0){

                    sCurrent = bestIteration;

                    fixedBlock[bestBlock] = 1;
                }
            }
        }

        end_timeMH = get_time_in_seconds();
        currentTime = end_timeMH - start_timeMH;

        UpdatePoolSolutions(bestPath, method, runData.debug);
    }
}

#endif
