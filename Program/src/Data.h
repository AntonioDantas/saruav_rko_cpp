#ifndef _DATA_H
#define _DATA_H

struct TSol
{
    std::vector <double> rk;
    std::vector<std::vector<int>> routes;
    double ofv = std::numeric_limits<double>::infinity();
    double f1;
    double f2;

    double best_time = 0.0;
    char nameMH[256];
};

struct TRunData
{
    int strategy;
    int control;
    int MAXTIME;
    int MAXRUNS;
    int debug;
    float restart;
    int sizePool;
    int VMIN;
    int VMAX;
    double F1MAX;
    double F2MAX;
};

struct TQ
{
    int S;
    double pVar;
    double q;
    int k;
    int kImp;
};

struct TState
{
    int label;
    std::vector<double> par;
    double ci;
    int numN;
    std::vector<int> Ai;
    std::vector<double> Qa;
    double sumQ;
    double maxQ;
    int minA;
    int maxA;
};

#endif
