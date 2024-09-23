/********************************************************
 * MC-MOSA method:
 * Multi-case Multi-Objective Simulated Annealing method 
 * ******************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "define.h"
#include "common1.h"
#include "common2.h"
#include "deviance.h"
#include "mcmosa.h"
#include "reaxffio.h"
#include "pso.h"

/*double RandomMove(){
	int i;
	for(i=0; i<nparas; i++){
		if ((double)rand()/RAND_MAX <= pcutoff)
			dp[i] = maxdp[i] * (2*((double)rand()/RAND_MAX) - 1);
		else
			dp[i] = 0.0;
		*para[i] += dp[i];	
		if( (*para[i]>uplimit[i])||(*para[i]<downlimit[i]) )
			{*para[i] -= dp[i]; dp[i]=0.0;}
	}
	return 0;
}
*/

double McmosaCompare(double *devi_new, double *dev_old){

	int i;
	double p; //possibility of the acceptance 
	double a, b, var, total_devi;

	p=0.0;
	total_devi=0.0;
	for(i=0; i<nobjs; i++){
		a = devi_old[i];
		b = devi_new[i];
		total_devi += lamda[i] * b;
	}
	var = exp(-1*total_devi/(coeff*temp));
	p = min(1, var);
	return p;
}


int McmosaUpdata(double p, double *devi){
	int i, flag, option=0;
	flag = 0;
	if(((double)rand()/RAND_MAX)<=p)
		flag = 1;

	/* Added new solution into paretos */
	if(flag){
		DeviCopy(devi, devi_old);	//Copy devi_new to devi_old
		UpdateParetos(devi);
	}
	else{
		//the movement is not accepted, delete random movements of paras
		for(i=0; i<nparas; i++)
			*para[i] -= dp[i];
	}

	//UpdatePSO();
	/* Find the best result in paretos every LOOP_INTERVAL stels */
	if((nloop%interval1)==0)
		option = 1;
	FindBestPareto(option);
	
	return flag;
}


int UpdateParetos(double *devi){
	int i;
	double var=0.0;

	for(i=0; i<nobjs; i++){
		paretos[nParetos].devi[i] = devi[i];
		var += lamda[i] * devi[i];
	}
	paretos[nParetos].total_devi = var;

	for(i=0; i<nparas; i++)	
		paretos[nParetos].para[i] = *para[i];
	nParetos++;
	
	while(nParetos>=MAX_PARETOS)
		DeleteBadPareto();
	return 0;
}

int DeleteBadPareto(){
	
	int i, j, remark=0;
	double var;

	var = 0.0;
	for(i=0; i<nParetos; i++)
		if(paretos[i].total_devi>var){
			var = paretos[i].total_devi;
			remark = i;
		}
	for(i=remark; i<nParetos; i++)
		if(i!=(nParetos-1)){
			for(j=0; j<nparas; j++)
				paretos[i].para[j] = paretos[i+1].para[j];
			for(j=0; j<nobjs;  j++)
				paretos[i].devi[j] = paretos[i+1].devi[j];
			paretos[i].total_devi = paretos[i+1].total_devi;
		}
	for(j=0; j<nparas; j++)
		paretos[nParetos-1].para[j] = 0.0;
	for(j=0; j<nobjs; j++)
		paretos[nParetos-1].devi[j] = 0.0;
	paretos[nParetos-1].total_devi = 0.0;
	nParetos--;

	return 0;
}


int FindBestPareto(int option){

	int    i, obj, remarkg, *remarkp;
	double varg, *varp;

	varg = 1.0e20;
	varp = (double *) malloc(sizeof(double)*nobjs);
	remarkg = 0;
	remarkp = (int *) malloc(sizeof(int)   *nobjs);
	for(obj=0; obj<nobjs; obj++){
		varp[obj] = 1.0e20;
		remarkp[obj] = 0;
	}

	for(i=0; i<nParetos; i++){
		if(paretos[i].total_devi<varg){
			varg = paretos[i].total_devi;
			remarkg = i;
		}
		for(obj=0; obj<nobjs; obj++){
			if(paretos[i].devi[obj] < varp[obj]){
				varp[obj] = paretos[i].devi[obj];
				remarkp[obj] = i;
			}
		}
	}

	for(i=0; i<nparas; i++){
		gbest[i] = paretos[remarkg].para[i];
		for(obj=0; obj<nobjs; obj++)
			pbest[obj][i] = paretos[remarkp[obj]].para[i];
	}
	if(option==1){
		for(i=0; i<nparas; i++)
			*para[i] = paretos[remarkg].para[i];
		for(i=0; i<nobjs; i++)
			devi_new[i] = paretos[remarkg].devi[i];	
		DeviCopy(devi_new, devi_old);	//Copy devi_new to devi_old
	}
	return 0;
}


int CoolTemp(){
	if(cool_mode==0){
		if(temp>MiddleTemp)
			temp *= coolp1;
		else
			temp -= coolp2;	
	}
	else if(cool_mode==1){
		if(temp>MiddleTemp)
			temp *= coolp1;
		else
			temp *= coolp2;
	}
	else if(cool_mode==2){
		if(temp>MiddleTemp)
			temp -= coolp1;
		else
			temp -= coolp2;
	}
	return 1;
}

