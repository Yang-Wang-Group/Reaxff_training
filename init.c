#define EXTERN extern
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <math.h>
#include "init.h"
#include "mcmosa.h"
#include "output.h"
#include "reaxffio.h"
#include "define.h"
#include "common1.h"
#include "common2.h"



int InitSys(){
	InitiationPartOne();
	Initiation4reaxff();
	ReadFiles();		
	InitiationPartTwo();
	InitiationPSO();
	calcu_Eta();
	ReadInputInfo();
	LamdaAssign();
	WriteOutput(0);
	InitCalc();
	return 1;
}


int InitiationPartOne(){

	int i, j, obj;

	/* MCMO-SA part */
	nparas = 0;
	nobjs  = 0;
	nParetos = 0;
	nloop  = 0;
	for(i=0; i<MAX_PARAMETERS; i++) {
		maxdp[i] = 0.0;
		dp[i] = 0.0;
	}
	for(obj=0; obj<MAX_TARGETS; obj++){
		opt[obj].npoints_ref = 0;
		opt[obj].npoints_md  = 0;
		for(i=0; i<MAX_DATAPOINTS; i++){
			opt[obj].ref_axisx[i] = 0.0;
			opt[obj].ref_axisy[i] = 0.0;
			opt[obj].atomicity[i] = 0;
			opt[obj].n_atom_species[i] = 0;
			for(j=0; j<MAX_ATOMICITY; j++){
				opt[obj].atom_type[i][j] = 0;
				opt[obj].posx[i][j] = 0.0;
				opt[obj].posy[i][j] = 0.0;
				opt[obj].posz[i][j] = 0.0;
			}
			for(j=0; j<3; j++)
				opt[obj].box[i][j] = 0.0;
		}
	}
	for(i=0; i<MAX_PARETOS; i++){
		for(j=0; j<MAX_PARAMETERS; j++)
			paretos[i].para[j] = 0.0;
		for(j=0; j<MAX_TARGETS; j++)
			paretos[i].devi[j] = 0.0;
		paretos[i].total_devi = 0.0;
		paretos[i].deleted_remark = 0;
	}

	/* Input Information */
	InitTemp = 8000.0;  MiddleTemp = 1000.0; FinalTemp = 5.0;
	cool_mode = 0; coolp1 = 0.999; coolp2 = 1.0;
	interval1 = 500; interval2 = 100;
	lamda_mode = 0; sprintf(lamdafile, "assign_lamda");
	pcutoff = 0.5; coeff = 1.0;
	sprintf(exename,"lmp_mpi");
	sprintf(outname,"log.lammps");
	return 0;
}


int ReadFiles(){
	//filename is "file-XXX-AB"
	//for example, file.bond.H-O, file.angle.H-O-C, file.EOS.C, file.EOS.O-Si

	int i, j, id;
	double a, b, c;
	char buff[512], na[2], fname[128], str1[30], str2[30];
	int flag, ip=0, obj; 
	DIR *dirp;
	FILE *fp;
	struct dirent *dp;

	dirp = opendir(".");
	obj = 0;
	opt[obj].min_energy_axisy = 20.0;
	while ( (dp=readdir(dirp))!=NULL ){
		sprintf(fname, "%s", dp->d_name);
		if (!strncmp(fname, "file", 4)){
			sscanf(fname, "%*[^-]-%[^-]-%s", str1, str2);
			sprintf(opt[obj].opt_type, "%s", str1);
			sprintf(opt[obj].opt_detail, "%s", str2);
			//opt_type indicates bond energy, angle energy, or others
			//opt_detail indicates which kinds of atoms containing in the system
			
			flag=0; ip=0;
			fp = fopen(fname,"r");
			for( ; (fgets(buff, 512, fp)!=NULL); ){
				if (!strncmp(buff, "Reference", 9)){
					sscanf(buff, "%*s %d %*s", &i);
					opt[obj].npoints_ref = i;
					for(i=0; i<opt[obj].npoints_ref; i++){
						fgets(buff, 512, fp);
						sscanf(buff, "%lf %lf %*s", &a, &b);
						opt[obj].ref_axisx[i] = a;
						opt[obj].ref_axisy[i] = b;
						if(b < opt[obj].min_energy_axisy){
							opt[obj].min_energy_axisx = a;
							opt[obj].min_energy_axisy = b;
						}
					}
				}
				if (!strncmp(buff, "Structure", 9)){
					sscanf(buff, "%*s %d %*s", &i);
					opt[obj].npoints_md = i;
					flag = 1; continue;}
				if (flag){
					if(!strncmp(buff,"Box", 3)){
						sscanf(buff, "%*s %lf %lf %lf", &a, &b, &c);
						opt[obj].box[ip][0] = a;
						opt[obj].box[ip][1] = b;
						opt[obj].box[ip][2] = c;
					}
					if(!strncmp(buff, "Atom", 4)){
						sscanf(buff, "%*s %d %d %*s", &i, &j); 
						opt[obj].atomicity[ip] = i;
						opt[obj].n_atom_species[ip] = num_ovar;
						for(id=0; id<opt[obj].atomicity[ip]; id++){
							fgets(buff, 512, fp);
							sscanf(buff, "%s %lf %lf %lf", na, &a, &b, &c);
							opt[obj].atom_type[ip][id] = find_atom_type(na);
							opt[obj].posx[ip][id] = a;
							opt[obj].posy[ip][id] = b;
							opt[obj].posz[ip][id] = c;
						}
						ip++;
					}
				}
			}
			if (opt[obj].npoints_md!=ip){
				printf("Error in %s: number of points are not same for MD calculation\n", fname); 
				printf("%d\t%d\n",opt[obj].npoints_md, ip); exit(0);
			}

			/* when opt object is bond, EOS, angle, or torsion, npoints_md must equal to npoints_ref */
			if ( (!strncmp(opt[obj].opt_type, "bond", 4)) || (!strncmp(opt[obj].opt_type, "EOS", 3)) || (!strncmp(opt[obj].opt_type, "eos", 3)) || (!strncmp(opt[obj].opt_type, "angle", 5)) || (!strncmp(opt[obj].opt_type, "torsion", 7)) )
				if (opt[obj].npoints_ref != opt[obj].npoints_md) 
					{printf("Error!!! npoints_ref != npoints_md for bond|EOS|angle|torsion optimization\n"); exit(0);}

			fclose(fp);
			obj++;
		}
	}
	nobjs = obj; //nobjs is the number of objectives in the error evaluation, equal to the number of file.XXX.A-B files.
	return 0; 
}

/*Calculate the attention weight for each objective*/
int calcu_Eta(){
	int obj,i,ip;
	double a,distance;
	double k[50];

	for (obj=0; obj<nobjs; obj++){
		a=0.0;
			if     (!strncmp(opt[obj].opt_type, "bond", 4)){
				for (ip=0; ip<opt[obj].npoints_ref; ip++){
				distance = fabs(opt[obj].ref_axisx[ip] - opt[obj].min_energy_axisx);
					if (distance < 2.0)
					k[ip] = 0.7-0.2*distance;
					else 
					k[ip] = 0.1;
				a += k[ip];
				}
				for(ip=0;ip<opt[obj].npoints_ref;ip++) opt[obj].eta[ip] = k[ip];
			}
			else if( (!strncmp(opt[obj].opt_type, "EOS", 3))||(!strncmp(opt[obj].opt_type, "eos", 3)) ){
				for(ip=0; ip<opt[obj].npoints_ref; ip++){
					distance = fabs(opt[obj].ref_axisx[i] - opt[obj].min_energy_axisx);
					if (distance < 1.0)
					k[ip] = 0.9-0.1*distance;
					else 
					{k[ip] = 0.5;}
					a += k[ip];
					}
				for (i=0; i<opt[obj].npoints_ref; i++)	opt[obj].eta[i] = k[i]/a; 
			}
			else if(!strncmp(opt[obj].opt_type, "angle", 5)){
				for(ip=0; ip<opt[obj].npoints_ref; ip++){
					distance = fabs(opt[obj].ref_axisx[ip] - opt[obj].min_energy_axisx);
					if (distance < 60)
					k[ip] = 0.25-distance/360;
					else 
					k[ip] = 0.083333;
				a += k[ip];
				}
				for(ip=0;ip<opt[obj].npoints_ref;ip++) opt[obj].eta[ip] = k[ip];
			}
	}
	return 0;
}

int ReadInputInfo(){

	int i=0, stopflag;
	double var;
	char buff[512];
	FILE *fp;
		
	fp = fopen("inputinfo","r");
	if(!fp) {printf("There is no inputinfo file!\n"); exit(0);}

	stopflag = 0;
	while(!stopflag){
		fgets(buff, 512, fp);
		if( !strncmp(buff, "temp", 4) )
			sscanf(buff, "%*s %lf %lf %lf", &InitTemp, &MiddleTemp, &FinalTemp);
		if( !strncmp(buff, "cool", 4) )
			sscanf(buff, "%*s %d %lf %lf", &cool_mode, &coolp1, &coolp2);
		if( !strncmp(buff, "maxloop", 7) )
			sscanf(buff, "%*s %d %d", &interval1, &interval2);
		if( !strncmp(buff, "lamda", 5) )
			sscanf(buff, "%*s %d %s", &lamda_mode, lamdafile);
		if( !strncmp(buff, "pcutoff", 7) )
			sscanf(buff, "%*s %lf", &pcutoff);
		if( !strncmp(buff, "coeff", 5) )
			sscanf(buff, "%*s %lf", &coeff);
		if( !strncmp(buff, "lammps", 6) )
			sscanf(buff, "%*s %s", exename);
		if( !strncmp(buff, "outname", 7) )
			sscanf(buff, "%*s %s", outname);
		if( !strncmp(buff, "pso", 3) ){
			sscanf(buff, "%*s %lf %lf %lf %lf", &comega, &cg, &var, &cnoise);
			for(i=0; i<nobjs; i++) cp[i] = var/(double)nobjs;
		}
		if( !strncmp(buff, "end", 3)||!strncmp(buff, "###", 3) )
			stopflag = 1;
		i++;
		if(i>1000) {printf("Error in ReadInputInfo()!!\n"); exit(0);}
	}
	sprintf(run, "%s > %s", exename, outname);
	return 1;		
}


double find_atom_weight(char type[2]){
	double a=0.0;
	if      (!strcmp(type, "X ")) {a=100.00000; return a;}

	/* Periodic 1, 2, & 3 */
	else if (!strcmp(type, "H"))  {a=  1.00794; return a;}
	else if (!strcmp(type, "He")) {a=  4.00260; return a;}
	else if (!strcmp(type, "Li")) {a=  6.94100; return a;}
	else if (!strcmp(type, "Be")) {a=  9.01218; return a;}
	else if (!strcmp(type, "B")) {a= 10.81100; return a;}
	else if (!strcmp(type, "C")) {a= 12.01100; return a;}
	else if (!strcmp(type, "N")) {a= 14.00674; return a;}
	else if (!strcmp(type, "O")) {a= 15.99940; return a;}
	else if (!strcmp(type, "F")) {a= 18.99840; return a;}
	else if (!strcmp(type, "Ne")) {a= 20.17970; return a;}
	else if (!strcmp(type, "Na")) {a= 22.989768;return a;}
	else if (!strcmp(type, "Mg")) {a= 24.30500; return a;}
	else if (!strcmp(type, "Al")) {a= 26.981539;return a;}
	else if (!strcmp(type, "Si")) {a= 28.08550; return a;}
	else if (!strcmp(type, "P")) {a= 30.973762;return a;}
	else if (!strcmp(type, "S")) {a= 32.06600; return a;}
	else if (!strcmp(type, "Cl")) {a= 35.45270; return a;}
	else if (!strcmp(type, "Ar")) {a= 39.94800; return a;}

	/* Periodic 4 */
	else if (!strcmp(type, "K")) {a= 39.09830; return a;}
	else if (!strcmp(type, "Ca")) {a= 40.07800; return a;}
	else if (!strcmp(type, "Sc")) {a= 44.95591; return a;}
	else if (!strcmp(type, "Ti")) {a= 47.86700; return a;}
	else if (!strcmp(type, "Ba")) {a= 50.94150; return a;}
	else if (!strcmp(type, "Cr")) {a= 51.99610; return a;}
	else if (!strcmp(type, "Mn")) {a= 54.93805; return a;}
	else if (!strcmp(type, "Fe")) {a= 55.84700; return a;}
	else if (!strcmp(type, "Co")) {a= 58.93320; return a;}
	else if (!strcmp(type, "Ni")) {a= 58.69340; return a;}
	else if (!strcmp(type, "Cu")) {a= 63.54600; return a;}
	else if (!strcmp(type, "Zn")) {a= 65.39000; return a;}
	else if (!strcmp(type, "Ga")) {a= 69.72300; return a;}
	else if (!strcmp(type, "Ge")) {a= 72.61000; return a;}
	else if (!strcmp(type, "As")) {a= 74.92159; return a;}
	else if (!strcmp(type, "Se")) {a= 78.96000; return a;}
	else if (!strcmp(type, "Br")) {a= 79.90400; return a;}
	else if (!strcmp(type, "Kr")) {a= 83.79800; return a;}

	/* Periodic 5 */
	else if (!strcmp(type, "Rb")) {a= 85.46800; return a;}
	else if (!strcmp(type, "Sr")) {a= 87.62000; return a;}
	else if (!strcmp(type, "Y")) {a= 88.90600; return a;}
	else if (!strcmp(type, "Zr")) {a= 91.22400; return a;}
	else if (!strcmp(type, "Nb")) {a= 92.90600; return a;}
	else if (!strcmp(type, "Mo")) {a= 95.95000; return a;}
	else if (!strcmp(type, "Tc")) {a= 98.00000; return a;}
	else if (!strcmp(type, "Ru")) {a=101.07000; return a;}
	else if (!strcmp(type, "Rh")) {a=102.91000; return a;}
	else if (!strcmp(type, "Pd")) {a=106.42000; return a;}
	else if (!strcmp(type, "Ag")) {a=107.87000; return a;}
	else if (!strcmp(type, "Cd")) {a=112.41000; return a;}
	else if (!strcmp(type, "In")) {a=114.82000; return a;}
	else if (!strcmp(type, "Sn")) {a=118.71000; return a;}
	else if (!strcmp(type, "Sb")) {a=121.76000; return a;}
	else if (!strcmp(type, "Te")) {a=127.60000; return a;}
	else if (!strcmp(type, "I ")) {a=126.90000; return a;}
	else if (!strcmp(type, "Xe")) {a=131.29000; return a;}

	/* Periodic 6 and 7 */
	else if (!strcmp(type, "W")) {a=183.84000; return a;}
	else if (!strcmp(type, "Pt")) {a=195.08000; return a;}
	else if (!strcmp(type, "Au")) {a=196.97000; return a;}
	else if (!strcmp(type, "Hg")) {a=200.59000; return a;}
	else if (!strcmp(type, "Pb")) {a=207.20000; return a;}
	else if (!strcmp(type, "Bi")) {a=208.98000; return a;}
	return a;
} 

char find_atom_type(char an[2]){
	int i, att=0;
	for(i=0; i<num_ovar; i++){
		if(strcmp(an, atomtype[i])==0)
			att = i + 1;
	}
	return att;	
}


int InitiationPartTwo(){
	int i;
	int mode_lamda;	//mod_lamda indicate how to initiate the lamda
			//lamda is the weight of each target in the optimization
	
	mode_lamda = 0;

	devi_new = (double *) malloc(sizeof(double)*nobjs);
	devi_old = (double *) malloc(sizeof(double)*nobjs);
	lamda    = (double *) malloc(sizeof(double)*nobjs);	
	origin_deviation = (double *) malloc(sizeof(double)*nobjs);
	for(i=0; i<nobjs; i++){
		devi_new[i] = 0.0;
		devi_old[i] = 0.0;
		lamda[i]    = 0.0;
		origin_deviation[i] = 0.0;
	}

	paraflag = (short *) malloc(sizeof(short)*nparas);
	for(i=0; i<nparas; i++){
		paraflag[i] = 0;
		if(maxdp[i]>0) paraflag[i] = 1;
	} 

	return 1;
}

int InitiationPSO(){
	int i, j;
	/* PSO initiation */
	comega = 0.0;
	cnoise = 0.0;
	cg    = 0.0;
	cp    = (double  *) malloc(sizeof(double  )*nobjs);
	pbest = (double **) malloc(sizeof(double *)*nobjs);	
	for(i=0; i<nobjs; i++){
		cp[i] = 0.0; //pso_cp/(double)nobjs;
		pbest[i] = (double *) malloc(sizeof(double)*nparas);
		for(j=0; j<nparas; j++) pbest[i][j] = 0.0;
	}
	gbest = (double *) malloc(sizeof(double)*nparas);
	for(i=0; i<nparas; i++) gbest[i] = 0.0;

	/* Initiation of dp (initial velocities in PSO) */
	for(i=0; i<nparas; i++) 
		dp[i] = maxdp[i] * (2*(double)rand()/RAND_MAX -1.0);

	return 1;
}


int LamdaAssign(){
	int i;
	if(lamda_mode==0){
		for(i=0;i<nobjs;i++){
			lamda[i]=1.0/nobjs;
		}
	}
	return 1;
}


int InitCalc(){
	double total_var;
	total_var = DeviCalc(devi_new,origin_deviation);
	McmosaUpdata(1.0, devi_new);
	return 1;
}

