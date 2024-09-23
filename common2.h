#include "define.h"

/* Used in MCMO-SA method */
extern int    nparas, nobjs, nParetos, nloop;
extern double temp;
extern double *para[MAX_PARAMETERS], maxdp[MAX_PARAMETERS], dp[MAX_PARAMETERS], uplimit[MAX_PARAMETERS], downlimit[MAX_PARAMETERS];
extern double *devi_new, *devi_old, *lamda, *origin_deviation;
//extern double p1[5], p2[5], p3[5];
extern OptTargets opt[MAX_TARGETS];
extern Pareto paretos[MAX_PARETOS];

/* Input Information */
extern int    cool_mode, lamda_mode, interval1, interval2;
extern double InitTemp, MiddleTemp, FinalTemp;
extern double coolp1, coolp2, pcutoff, coeff;
extern char   lamdafile[20], exename[20], outname[20], run[40];

/* ReaxFF parameter part */
extern double gvar[N_GVAR];
extern double ovar[MAX_ATOMTYPES][N_OVAR];
extern double pvar[100][N_PVAR];
extern double vdwvar[100][N_VDWVAR];
extern double tvar[200][N_TVAR];
extern double fvar[100][N_FVAR];
extern double hvar[20][N_HVAR];
extern int    num_gvar, num_ovar, num_pvar, num_vdwvar, num_tvar, num_fvar, num_hvar;
extern char   atomtype[MAX_ATOMTYPES][3];
extern double atomweight[MAX_ATOMTYPES]; 
extern char   pname[100][20], vdwname[100][20], tname[200][20], fname[100][20], hname[20][20];
extern short  *paraflag;

/* PSO */
extern double **pbest, *gbest, *cp, cnoise, cg, comega;

