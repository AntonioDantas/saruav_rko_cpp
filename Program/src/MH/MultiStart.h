#ifndef _MultiStart_H
#define _MultiStart_H

void MultiStart(const TRunData &runData, const TProblemData &data)
{
    const char* method = "MultiStart";

    static TSol s;
    static TSol sBest;

    static int IterT;

    float currentTime = 0;
    IterT = 0;

    double start_timeMH = get_time_in_seconds();
    double end_timeMH = get_time_in_seconds();

    CreateInitialSolutions(s, data.n);
    s.ofv = Decoder(s, data);
    sBest = s;

    while(currentTime < runData.MAXTIME*runData.restart)
    {
        if (stop_execution.load()) return;

		CreateInitialSolutions(s, data.n);
        s.ofv = Decoder(s, data);

		if (s.ofv < sBest.ofv)
		{
			sBest = s;

            UpdatePoolSolutions(sBest, method, runData.debug);
		}

        IterT++;

        end_timeMH = get_time_in_seconds();
        currentTime = end_timeMH - start_timeMH;
    }
}

#endif
