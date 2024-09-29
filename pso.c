/********************************************************
 * Particle Swarm Optimization (PSO):
 * ******************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "define.h"
#include "common1.h"
#include "common2.h"
#include "pso.h"

double RandomMove(){
	int i;
	srand((unsigned)time(NULL));
	for(i=0; i<nparas; i++){
		varflag[i] = 0;
		if(paraflag[i]&1){
			if ((double)rand()/RAND_MAX <= pcutoff)
			dp[i] = PSO(i);
			*para[i] += dp[i];
			varflag[i] = 1;
			if ( (*para[i]>uplimit[i]) || (*para[i]<downlimit[i])){
				*para[i] -= dp[i]; 
				varflag[i] = 0;
			}

		}
	}
	return 0;
}

double PSO(int i){
	int obj;
	double val = 0;	

	val += comega * dp[i];
	for (obj=0; obj<nobjs; obj++)
		val += cp[obj] * (pbest[obj][i] - *para[i]) * (double)rand()/RAND_MAX; 

	val += cg * (gbest[i] - *para[i]) * (double)rand()/RAND_MAX;
	val += cnoise * maxdp[i] * (2*(double)rand()/RAND_MAX - 1.0);

	if (val > maxdp[i]) val = maxdp[i];
	return val;
}
