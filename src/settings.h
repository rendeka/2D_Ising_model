#pragma once
#define N_SWEEP 500 //number of sweeps
#define N_TSWEEP 50 //number of sweeps for thermalization
#define N_INDEP 100
#define n 16 //length of the grid
#define N_THREAD 8 //number of threads
#define J 1 //magnetic bonding constant
#define nT 30 //number of different temperatures for simalation
#define Tfirst 1.5 //lower temperature we simulate
#define Tlast 3.0 //highest temperature we simulate
#define Tstep (Tlast - Tfirst) / nT
#define normalization double(N_SWEEP) * double(n) * double(n) * double(N_INDEP) //normalization constant
