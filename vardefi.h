#ifndef _VARDEFI_H_
#define _VARDEFI_H_

#include "define.h"

/* Used in MCMO-SA method */
int    nparas, nobjs, nParetos, nloop;
double temp;
double *para[MAX_PARAMETERS], maxdp[MAX_PARAMETERS], dp[MAX_PARAMETERS], uplimit[MAX_PARAMETERS], downlimit[MAX_PARAMETERS];
double *devi_new, *devi_old, *lamda, *origin_deviation;
OptTargets opt[MAX_TARGETS];
Pareto paretos[MAX_PARETOS];

/* Input Information */
int    cool_mode, lamda_mode, interval1, interval2;
double InitTemp, MiddleTemp, FinalTemp;
double coolp1, coolp2, pcutoff, coeff;
char   lamdafile[20], exename[20], outname[20], run[40];

/* ReaxFF parameter part */
double gvar[N_GVAR];		//global parameters
double ovar[MAX_ATOMTYPES][N_OVAR];	//one-atom parameters
double pvar[100][N_PVAR];	//pair-atom parameters
double vdwvar[100][N_VDWVAR];	//vdw parameters
double tvar[200][N_TVAR];	//triple-atom parameters
double fvar[100][N_FVAR];	//four-atom parameters
double hvar[20][N_HVAR];	//hydrogen bond parameters
int    num_gvar, num_ovar, num_pvar, num_vdwvar, num_tvar, num_fvar, num_hvar;
char   atomtype[MAX_ATOMTYPES][3];
double atomweight[MAX_ATOMTYPES]; 
char   pname[100][20], vdwname[100][20], tname[200][20], fname[100][20], hname[20][20];
short  *paraflag;

/* PSO */
double **pbest, *gbest, *cp, cnoise, cg, comega;

#endif
