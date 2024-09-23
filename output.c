#include <stdio.h>
#include "output.h"
#include "define.h"
#include "common1.h"
#include "common2.h"
#include "reaxffio.h"


int WriteOutput(int mode){
//If mode==0, open output file and write the head line
//If mode>0(i.e mode==1), write output lines every step
	WriteParameters(mode);
	WriteDeviations(mode);
	WriteOriginDeviations(mode);
	return 0;
}


int WriteParameters(int mode){

	FILE *fp;

	if(mode==0){
		fp = fopen("result_parameters","w");
		fprintf(fp, "temp\tparameters....\n");
		fclose(fp);
	}
	else if(mode){
		fp = fopen("result_parameters", "a");
		write_para_everystep(fp);
		fclose(fp);
	}
	return 0;
}


int WriteDeviations(int mode){
	int i, obj;
	double total_devi=0.0;
	FILE *fp;

	if(mode==0){
		fp = fopen("result_deviance","w");
		fprintf(fp, "temp\t");
		for(obj=0; obj<nobjs; obj++)
			fprintf(fp, "%s-%s\t", opt[obj].opt_type, opt[obj].opt_detail); 
		fprintf(fp,"total_deviance\n");
		fclose(fp);
		WriteEta();
	}
	else if(mode){
		fp = fopen("result_deviance", "a");
		fprintf(fp, "%lf", temp);
		for(i=0; i<nobjs; i++){
			fprintf(fp, "\t%lf",devi_new[i]);
			total_devi += lamda[i] * devi_new[i];
		}		
		fprintf(fp, "\t%lf\n",total_devi);
		fclose(fp);
	}
	return 0;
}

int WriteOriginDeviations(int mode){
	int i, obj;
	double total_devi=0.0;
	FILE *fp;

	if(mode==0){
		fp = fopen("result_origin_deviance","w");
		fprintf(fp, "temp\t");
		for(obj=0; obj<nobjs; obj++)
			fprintf(fp, "%s-%s\t", opt[obj].opt_type, opt[obj].opt_detail); 
		fprintf(fp,"total_origin_deviance\n");
		fclose(fp);
	}
	else if(mode){
		fp = fopen("result_origin_deviance", "a");
		fprintf(fp, "%lf", temp);
		for(i=0; i<nobjs; i++){
			fprintf(fp, "\t%lf",origin_deviation[i]);
			total_devi += lamda[i] * origin_deviation[i];
		}		
		fprintf(fp, "\t%lf\n",total_devi);
		fclose(fp);
	}
	return 0;
}

int WriteEta(){
	int i,obj;
	FILE *fp;
	fp= fopen("Eta.txt","w");
	for(obj=0;obj<nobjs;obj++){
		fprintf(fp,"%d\n",obj);
		for(i=0;i<opt[obj].npoints_ref;i++){
			fprintf(fp,"\n%lf\t",opt[obj].eta[i]);
		}
	}
	fclose(fp);
	return 0;
}

int WriteOutputParetos(){
	int i, j;
	FILE *fp;
	fp = fopen("result_paretos","w");

	//head line
	fprintf(fp, "#:\t%d devis, total_devi, %d paras\n",nobjs,nparas);

	//main lines of result_paretos	
	for(i=0; i<nParetos; i++){
		for(j=0; j<nobjs; j++)
			fprintf(fp, "\t%lf", paretos[i].devi[j]);
		fprintf(fp, "%d\t%lf", i, paretos[i].total_devi);
		for(j=0; j<nparas; j++)
			fprintf(fp, "\t%lf", paretos[i].para[j]);
		fprintf(fp, "\n");
	}
	fclose(fp);

	return 1;
}


