#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "common1.h"
#include "vardefi.h" 
#include "common2.h"
#include "reaxffio.h"
#include "init.h"

double  Initiation4reaxff(){
	//read reaxff parameter file
	ReadParas4reaxff();
	//read max displacement for each parameter
	ReadMaxdp4reaxff();
	//read uplimit value for each parameter file
	ReadUplimit4reaxff();
	//read downlimit value for each parameter file
	ReadDownlimit4reaxff();
	return 1;
}


double write_reaxff_input_files(int obj, int ip){
	//write_def_rd();
	write_data_md(obj, ip); 
	write_para_reax();       
	return 1;
}

/*
double write_input_rd(int obj, int ip){
	char linedata[200]={0};
	int i;
	FILE *fp;
	FILE *fpw;
	fp = fopen("input.in","r");
	fpw=fopen("input2.in","w");
	while (fgets(linedata,sizeof(linedata)-1,fp))
    {
        if (strcmp(linedata,"pair_coeff")==0)
        {
			fprintf(fpw,"%s %c %c ");
			for(i=0;i<opt[obj].atomicity[ip];i++){
            	fprintf(fpw,"%s ",pair_coeff,*,*,atomtype[i]);
			}
        }else
            fputs(linedata,fpw);            
    }
    fclose(fp);
    fclose(fpw);
    system("del data.txt");
    system("rename tmp.txt data.txt");
		
}*/

double write_para_reax(){

	int i;	
	FILE *fp;
	fp = fopen("para.reax", "w");

	fprintf(fp, "Reactive MD-force field: used for parameter fitting\n");
	/* Write global parameter part */
	fprintf(fp, " %d\t! Number of general parameters\n", num_gvar);
	for(i=0; i<num_gvar; i++)
		fprintf(fp, "%10.4f !XXX\n", gvar[i]);

	/* Write one-atom parameter part */
	fprintf(fp, "%3d   !Nr of atoms\n", num_ovar); 
	fprintf(fp, "       XXXXXXXXXXXXXXXX\n");
	fprintf(fp, "       XXXXXXXXXXXXXXXX\n");
	fprintf(fp, "       XXXXXXXXXXXXXXXX\n");
	for(i=0; i<num_ovar; i++){
		fprintf(fp,"%s   %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", atomtype[i], ovar[i][0], ovar[i][1], ovar[i][2], ovar[i][3], ovar[i][4], ovar[i][5], ovar[i][6], ovar[i][7]);
		fprintf(fp, "    %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",              ovar[i][8], ovar[i][9], ovar[i][10],ovar[i][11],ovar[i][12],ovar[i][13],ovar[i][14],ovar[i][15]);
		fprintf(fp, "    %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",              ovar[i][16],ovar[i][17],ovar[i][18],ovar[i][19],ovar[i][20],ovar[i][21],ovar[i][22],ovar[i][23]);
		fprintf(fp, "    %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",              ovar[i][24],ovar[i][25],ovar[i][26],ovar[i][27],ovar[i][28],ovar[i][29],ovar[i][30],ovar[i][31]);
	}	

	/* Write pair-atom parameter part */
	fprintf(fp, " %d      !Nr of bonds; XXX\n                     XXXXXXXXXXX\n", num_pvar);
	for(i=0; i<num_pvar; i++){
		fprintf(fp,     "%s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", pname[i], pvar[i][0], pvar[i][1], pvar[i][2], pvar[i][3], pvar[i][4], pvar[i][5], pvar[i][6], pvar[i][7]);
		fprintf(fp, "       %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",           pvar[i][8], pvar[i][9], pvar[i][10],pvar[i][11],pvar[i][12],pvar[i][13],pvar[i][14],pvar[i][15]);
	}

	/* Write vdw parameter part */
	fprintf(fp, "%3d    !Nr of off-diagonal terms;XXXXXX\n", num_vdwvar);
	for(i=0; i<num_vdwvar; i++)
		fprintf(fp, "%s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", vdwname[i], vdwvar[i][0], vdwvar[i][1], vdwvar[i][2], vdwvar[i][3], vdwvar[i][4], vdwvar[i][5]);

	/* Write triple-atom parameter part */
	fprintf(fp, "%3d    !Nr of angles;XXXXXX\n", num_tvar);
	for(i=0; i<num_tvar; i++)
		fprintf(fp, "%s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", tname[i], tvar[i][0], tvar[i][1], tvar[i][2], tvar[i][3], tvar[i][4], tvar[i][5], tvar[i][6]);

	/* Write four-atom parameter part */
	fprintf(fp, "%3d    !Nr of torsions; XXX\n", num_fvar);	
	for(i=0; i<num_fvar; i++)
		fprintf(fp, "%s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", fname[i], fvar[i][0], fvar[i][1], fvar[i][2], fvar[i][3], fvar[i][4], fvar[i][5], fvar[i][6]);

	/* Write hydrogen bond parameter part */
	fprintf(fp, "%3d    !Nr of hydrogen bonds; XXX\n", num_hvar);
	for(i=0; i<num_hvar; i++)
		fprintf(fp, "%s %8.4f %8.4f %8.4f %8.4f\n", hname[i], hvar[i][0], hvar[i][1], hvar[i][2], hvar[i][3]);
	fclose(fp);
	return 1;
}


double write_para_everystep(FILE *fp){

	int i;	

	/* Write head line */
	fprintf(fp, "Temp: %f --------------------------------------------\n", temp);

	/* Write global parameter part */
	fprintf(fp, " %d\t! Number of general parameters\n", num_gvar);
	for(i=0; i<num_gvar; i++)
		fprintf(fp, "%10.4f !%d\n", gvar[i], i+1);

	/* Write one-atom parameter part */
	fprintf(fp, "%3d  !Nr of one atom parameters\n", num_ovar); 
	for(i=0; i<num_ovar; i++){
		fprintf(fp,"%s   %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", atomtype[i], ovar[i][0], ovar[i][1], ovar[i][2], ovar[i][3], ovar[i][4], ovar[i][5], ovar[i][6], ovar[i][7]);
		fprintf(fp, "    %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",              ovar[i][8], ovar[i][9], ovar[i][10],ovar[i][11],ovar[i][12],ovar[i][13],ovar[i][14],ovar[i][15]);
		fprintf(fp, "    %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",              ovar[i][16],ovar[i][17],ovar[i][18],ovar[i][19],ovar[i][20],ovar[i][21],ovar[i][22],ovar[i][23]);
		fprintf(fp, "    %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",              ovar[i][24],ovar[i][25],ovar[i][26],ovar[i][27],ovar[i][28],ovar[i][29],ovar[i][30],ovar[i][31]);
	}	

	/* Write pair-atom parameter part */
	fprintf(fp, " %d      !Nr of pari atom parameters\n", num_pvar);
	for(i=0; i<num_pvar; i++){
		fprintf(fp,     "%s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", pname[i], pvar[i][0], pvar[i][1], pvar[i][2], pvar[i][3], pvar[i][4], pvar[i][5], pvar[i][6], pvar[i][7]);
		fprintf(fp, "       %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",           pvar[i][8], pvar[i][9], pvar[i][10],pvar[i][11],pvar[i][12],pvar[i][13],pvar[i][14],pvar[i][15]);
	}

	/* Write vdw parameter part */
	fprintf(fp, "%3d    !Nr of off-diagonal terms\n", num_vdwvar);
	for(i=0; i<num_vdwvar; i++)
		fprintf(fp, "%s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", vdwname[i], vdwvar[i][0], vdwvar[i][1], vdwvar[i][2], vdwvar[i][3], vdwvar[i][4], vdwvar[i][5]);

	/* Write triple-atom parameter part */
	fprintf(fp, "%3d    !Nr of 3-body parameters\n", num_tvar);
	for(i=0; i<num_tvar; i++)
		fprintf(fp, "%s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", tname[i], tvar[i][0], tvar[i][1], tvar[i][2], tvar[i][3], tvar[i][4], tvar[i][5], tvar[i][6]);

	/* Write four-atom parameter part */
	fprintf(fp, "%3d    !Nr of 4-body parameters\n", num_fvar);	
	for(i=0; i<num_fvar; i++)
		fprintf(fp, "%s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", fname[i], fvar[i][0], fvar[i][1], fvar[i][2], fvar[i][3], fvar[i][4], fvar[i][5], fvar[i][6]);

	/* Write hydrogen bond parameter part */
	fprintf(fp, "%3d    !Nr of hydrogen bonds; XXX\n", num_hvar);
	for(i=0; i<num_hvar; i++)
		fprintf(fp, "%s %8.4f %8.4f %8.4f %8.4f\n", hname[i], hvar[i][0], hvar[i][1], hvar[i][2], hvar[i][3]);

	fprintf(fp, "\n");
	return 1;
}


double analysis_for_reaxff_output(char filename[]){
	FILE *fp;
	int i;
	double a=0;
	char buff[512],*mark, *token;
	fp = fopen(filename, "r");
	if(fp){
		for(i=0; fgets(buff, 512, fp)!=NULL; i++){
			mark=strtok(buff," ");
			if(strcmp(mark,"Step")==0){ 
				fgets(buff, 512, fp); 
				token = strtok(buff, " ");
				token = strtok(NULL, " ");
				token = strtok(NULL, " ");	
				a = atof(token);
				break;
			}
		}
		fclose(fp);	
	}
	else{
		printf("No reaxff output file!\n");
		exit(1);
	}
	return a;
}


double ReadParas4reaxff(){

	int i, j, ipara, num;
	int at1, at2, at3, at4;
	double a, b, c, d, e, f, g, h;
	char atype[2];
	char buff[512];
	FILE *fp;	

	ipara=0;

	fp = fopen("para.reax", "r");
	if(!fp) exit(0);
	system("cp para.reax initpara.reax");
	
	//Read the first line
	fgets(buff, 512, fp);

	//Read the global parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);// Number of general parameters
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%lf %*s", &gvar[i]);
		//Following is used for the MC-MOSA calculation
		para[ipara] = &gvar[i];
		//maxdp[ipara]= 0.1*fabs(gvar[i]);
		ipara++;
	}	
	num_gvar = num;

	//Read the one-atom parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);//Nr of atoms
	fgets(buff, 512, fp); fgets(buff, 512, fp); fgets(buff, 512, fp);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%s %lf %lf %lf %lf %lf %lf %lf %lf", atype, &a, &b, &c, &d, &e, &f, &g, &h);
		atomtype[i][0] = atype[0];
		atomtype[i][1] = atype[1];
		atomweight[i]  = find_atom_weight(atomtype[i]);	
		ovar[i][0] = a; para[ipara] = &ovar[i][0]; ipara++;
		ovar[i][1] = b; para[ipara] = &ovar[i][1]; ipara++;
		ovar[i][2] = c; para[ipara] = &ovar[i][2]; ipara++;
		ovar[i][3] = d; para[ipara] = &ovar[i][3]; ipara++;
		ovar[i][4] = e; para[ipara] = &ovar[i][4]; ipara++;
		ovar[i][5] = f; para[ipara] = &ovar[i][5]; ipara++;
		ovar[i][6] = g; para[ipara] = &ovar[i][6]; ipara++;
		ovar[i][7] = h; para[ipara] = &ovar[i][7]; ipara++;
		for(j=1; j<=3; j++){
			fgets(buff, 512, fp);
			sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f, &g, &h);
			ovar[i][0+j*8] = a; para[ipara] = &ovar[i][0+j*8]; ipara++;
			ovar[i][1+j*8] = b; para[ipara] = &ovar[i][1+j*8]; ipara++;
			ovar[i][2+j*8] = c; para[ipara] = &ovar[i][2+j*8]; ipara++;
			ovar[i][3+j*8] = d; para[ipara] = &ovar[i][3+j*8]; ipara++;
			ovar[i][4+j*8] = e; para[ipara] = &ovar[i][4+j*8]; ipara++;
			ovar[i][5+j*8] = f; para[ipara] = &ovar[i][5+j*8]; ipara++;
			ovar[i][6+j*8] = g; para[ipara] = &ovar[i][6+j*8]; ipara++;
			ovar[i][7+j*8] = h; para[ipara] = &ovar[i][7+j*8]; ipara++;
		}
	}
	num_ovar = num;		

	//Read the pair-atom parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	fgets(buff, 512, fp); 
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &a, &b, &c, &d, &e, &f, &g, &h);
		sprintf(pname[i], " %2d %2d",at1,at2); 
		pvar[i][0] = a; para[ipara]=&pvar[i][0]; ipara++;
		pvar[i][1] = b; para[ipara]=&pvar[i][1]; ipara++;
		pvar[i][2] = c; para[ipara]=&pvar[i][2]; ipara++;
		pvar[i][3] = d; para[ipara]=&pvar[i][3]; ipara++;
		pvar[i][4] = e; para[ipara]=&pvar[i][4]; ipara++;
		pvar[i][5] = f; para[ipara]=&pvar[i][5]; ipara++;
		pvar[i][6] = g; para[ipara]=&pvar[i][6]; ipara++;
		pvar[i][7] = h; para[ipara]=&pvar[i][7]; ipara++;
		fgets(buff, 512, fp);
		sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f, &g, &h);
		pvar[i][8] = a; para[ipara]=&pvar[i][8]; ipara++;
		pvar[i][9] = b; para[ipara]=&pvar[i][9]; ipara++;
		pvar[i][10]= c; para[ipara]=&pvar[i][10];ipara++;
		pvar[i][11]= d; para[ipara]=&pvar[i][11];ipara++;
		pvar[i][12]= e; para[ipara]=&pvar[i][12];ipara++;
		pvar[i][13]= f; para[ipara]=&pvar[i][13];ipara++;
		pvar[i][14]= g; para[ipara]=&pvar[i][14];ipara++;
		pvar[i][15]= h; para[ipara]=&pvar[i][15];ipara++;
	}	
	num_pvar = num;
		
	//Read the vdw parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %lf %lf %lf %lf %lf %lf", &at1, &at2, &a, &b, &c, &d, &e, &f);
		sprintf(vdwname[i], " %2d %2d",at1,at2); 
		vdwvar[i][0] = a; para[ipara]=&vdwvar[i][0]; ipara++;
		vdwvar[i][1] = b; para[ipara]=&vdwvar[i][1]; ipara++;
		vdwvar[i][2] = c; para[ipara]=&vdwvar[i][2]; ipara++;
		vdwvar[i][3] = d; para[ipara]=&vdwvar[i][3]; ipara++;
		vdwvar[i][4] = e; para[ipara]=&vdwvar[i][4]; ipara++;
		vdwvar[i][5] = f; para[ipara]=&vdwvar[i][5]; ipara++;
	}
	num_vdwvar = num;
		
	//Read the triple-atom (angle) parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &at3, &a, &b, &c, &d, &e, &f, &g);
		sprintf(tname[i], " %2d %2d %2d",at1,at2,at3); 
		tvar[i][0] = a; para[ipara]=&tvar[i][0]; ipara++;
		tvar[i][1] = b; para[ipara]=&tvar[i][1]; ipara++;
		tvar[i][2] = c; para[ipara]=&tvar[i][2]; ipara++;
		tvar[i][3] = d; para[ipara]=&tvar[i][3]; ipara++;
		tvar[i][4] = e; para[ipara]=&tvar[i][4]; ipara++;
		tvar[i][5] = f; para[ipara]=&tvar[i][5]; ipara++;
		tvar[i][6] = g; para[ipara]=&tvar[i][6]; ipara++;
	}	
	num_tvar = num;

	//Read the four-atoms (torsion) parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &at3, &at4, &a, &b, &c, &d, &e, &f, &g);
		sprintf(fname[i], " %2d %2d %2d %2d",at1,at2,at3,at4); 
		fvar[i][0] = a; para[ipara]=&fvar[i][0]; ipara++;
		fvar[i][1] = b; para[ipara]=&fvar[i][1]; ipara++;
		fvar[i][2] = c; para[ipara]=&fvar[i][2]; ipara++;
		fvar[i][3] = d; para[ipara]=&fvar[i][3]; ipara++;
		fvar[i][4] = e; para[ipara]=&fvar[i][4]; ipara++;
		fvar[i][5] = f; para[ipara]=&fvar[i][5]; ipara++;
		fvar[i][6] = g; para[ipara]=&fvar[i][6]; ipara++;
	}	
	num_fvar = num;
		
	//Read the hydrogen bond parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %lf %lf %lf %lf", &at1, &at2, &at3, &a, &b, &c, &d);
		sprintf(hname[i], " %2d %2d %2d",at1,at2,at3); 
		hvar[i][0] = a; para[ipara]=&hvar[i][0]; ipara++;
		hvar[i][1] = b; para[ipara]=&hvar[i][1]; ipara++;
		hvar[i][2] = c; para[ipara]=&hvar[i][2]; ipara++;
		hvar[i][3] = d; para[ipara]=&hvar[i][3]; ipara++;
	}	
	num_hvar = num;

	fclose(fp);
	//Total number of parameters needed to be determine
	nparas = ipara;	

	return 0;
}


double ReadMaxdp4reaxff(){

	int i, j, ipara, num;
	int at1, at2, at3, at4;
	double a, b, c, d, e, f, g, h;
	char atype[2], name[20];
	char buff[512];
	FILE *fp;	

	ipara=0;

	fp = fopen("maxdp.reax", "r");
	if(!fp){
		printf("There is no maxdp.reax file!!!\n");
		exit(0);
	}
	
	//Read the first line
	fgets(buff, 512, fp);
		
	//Read the global parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%lf %*s", &a);
		maxdp[ipara] = fabs(a); ipara++;
	}
	if(num!=num_gvar) { printf("num_gvar is unconsistent! maxdp\n"); exit(0);}


	//Read the one-atom parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%s %lf %lf %lf %lf %lf %lf %lf %lf", atype, &a, &b, &c, &d, &e, &f, &g, &h);
		if( strcmp(atype, atomtype[i])!=0 ) {
			printf("%s   %s\n", atype, atomtype[i]);
			printf("atom type [%d] is unconsistent! maxdp\n",i); exit(0);}
		maxdp[ipara] = fabs(a); ipara++;
		maxdp[ipara] = fabs(b); ipara++;
		maxdp[ipara] = fabs(c); ipara++;
		maxdp[ipara] = fabs(d); ipara++;
		maxdp[ipara] = fabs(e); ipara++;
		maxdp[ipara] = fabs(f); ipara++;
		maxdp[ipara] = fabs(g); ipara++;
		maxdp[ipara] = fabs(h); ipara++;
		for(j=1; j<=3; j++){
			fgets(buff, 512, fp);
			sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f, &g, &h);
			maxdp[ipara] = fabs(a); ipara++;
			maxdp[ipara] = fabs(b); ipara++;
			maxdp[ipara] = fabs(c); ipara++;
			maxdp[ipara] = fabs(d); ipara++;
			maxdp[ipara] = fabs(e); ipara++;
			maxdp[ipara] = fabs(f); ipara++;
			maxdp[ipara] = fabs(g); ipara++;
			maxdp[ipara] = fabs(h); ipara++;
		}
	}
	if(num!=num_ovar) {printf("num_ovar is unconsistent! maxdp\n"); exit(0);}

	//Read the pair-atom parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num); 
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &a, &b, &c, &d, &e, &f, &g, &h);
		sprintf(name, " %2d %2d",at1,at2); 
		if( strcmp(name, pname[i])!=0 )
			{printf("pair-atom parameters of (%d,%d) are unconsistent! maxdp\n",at1,at2); exit(0);}
		maxdp[ipara]=fabs(a); ipara++;
		maxdp[ipara]=fabs(b); ipara++;
		maxdp[ipara]=fabs(c); ipara++;
		maxdp[ipara]=fabs(d); ipara++;
		maxdp[ipara]=fabs(e); ipara++;
		maxdp[ipara]=fabs(f); ipara++;
		maxdp[ipara]=fabs(g); ipara++;
		maxdp[ipara]=fabs(h); ipara++;
		fgets(buff, 512, fp);
		sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f, &g, &h);
		maxdp[ipara]=fabs(a); ipara++;
		maxdp[ipara]=fabs(b); ipara++;
		maxdp[ipara]=fabs(c); ipara++;
		maxdp[ipara]=fabs(d); ipara++;
		maxdp[ipara]=fabs(e); ipara++;
		maxdp[ipara]=fabs(f); ipara++;
		maxdp[ipara]=fabs(g); ipara++;
		maxdp[ipara]=fabs(h); ipara++;
	}	
	if(num!=num_pvar) {printf("num_pvar is unconsistent! maxdp\n"); exit(0);}

	//Read the vdw parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %lf %lf %lf %lf %lf %lf", &at1, &at2, &a, &b, &c, &d, &e, &f);
		sprintf(name, " %2d %2d",at1,at2); 
		if( strcmp(name, vdwname[i])!=0 )
			{printf("vdW parameters of (%d,%d) are unconsistent! maxdp\n",at1,at2); exit(0);}
		maxdp[ipara]=fabs(a); ipara++;
		maxdp[ipara]=fabs(b); ipara++;
		maxdp[ipara]=fabs(c); ipara++;
		maxdp[ipara]=fabs(d); ipara++;
		maxdp[ipara]=fabs(e); ipara++;
		maxdp[ipara]=fabs(f); ipara++;
	}
	if(num!=num_vdwvar) {printf("num_vdwvar is unconsistent! maxdp\n"); exit(0);}

	//Read the triple-atom (angle) parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &at3, &a, &b, &c, &d, &e, &f, &g);
		sprintf(name, " %2d %2d %2d", at1, at2, at3); 
		if( strcmp(name, tname[i])!=0 )
			{printf("triple-atom parameters of (%d,%d,%d) are unconsistent! maxdp\n",at1,at2,at3); exit(0);}
		maxdp[ipara]=fabs(a); ipara++;
		maxdp[ipara]=fabs(b); ipara++;
		maxdp[ipara]=fabs(c); ipara++;
		maxdp[ipara]=fabs(d); ipara++;
		maxdp[ipara]=fabs(e); ipara++;
		maxdp[ipara]=fabs(f); ipara++;
		maxdp[ipara]=fabs(g); ipara++;
	}	
	if(num!=num_tvar) {printf("num_tvar is unconsistent! maxdp\n"); exit(0);}
		
	//Read the four-atoms (torsion) parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &at3, &at4, &a, &b, &c, &d, &e, &f, &g);
		sprintf(name, " %2d %2d %2d %2d",at1,at2,at3,at4); 
		if( strcmp(name, fname[i])!=0 )
			{printf("4body parameters of (%d,%d,%d,%d) are unconsistent! maxdp\n",at1,at2,at3,at4); exit(0);}
		maxdp[ipara]=fabs(a); ipara++;
		maxdp[ipara]=fabs(b); ipara++;
		maxdp[ipara]=fabs(c); ipara++;
		maxdp[ipara]=fabs(d); ipara++;
		maxdp[ipara]=fabs(e); ipara++;
		maxdp[ipara]=fabs(f); ipara++;
		maxdp[ipara]=fabs(g); ipara++;
	}	
	if(num!=num_fvar) {printf("num_fvar is unconsistent! maxdp\n"); exit(0);}
		
	//Read the hydrogen bond parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %lf %lf %lf %lf", &at1, &at2, &at3, &a, &b, &c, &d);
		sprintf(name, " %2d %2d %2d",at1,at2,at3); 
		if( strcmp(name, hname[i])!=0 )
			{printf("hydrogen bond of (%d,%d,%d) are unconsistent! maxdp\n",at1,at2,at3); exit(0);}
		maxdp[ipara]=fabs(a); ipara++;
		maxdp[ipara]=fabs(b); ipara++;
		maxdp[ipara]=fabs(c); ipara++;
		maxdp[ipara]=fabs(d); ipara++;
	}	
	if(num!=num_hvar) {printf("num_hvar is unconsistent! maxdp\n"); exit(0);}
	
	fclose(fp);

	if(ipara!=nparas) {printf("Number of total parameters is not consistent to maxdp!\n"); exit(0);}

	return 0;
}


double ReadUplimit4reaxff(){

	int i, j, ipara, num;
	int at1, at2, at3, at4;
	double a, b, c, d, e, f, g, h;
	char atype[2], name[20];
	char buff[512];
	FILE *fp;	

	ipara=0;

	fp = fopen("uplimit.reax", "r");
	if(!fp){
		printf("There is no uplimit.reax file!!!\n");
		exit(0);
	}
	
	//Read the first line
	fgets(buff, 512, fp);
		
	//Read the global parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%lf %*s", &a);
		uplimit[ipara] = fabs(a); ipara++;
	}
	if(num!=num_gvar) { printf("num_gvar is unconsistent! uplimit\n"); exit(0);}

	//Read the one-atom parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%s %lf %lf %lf %lf %lf %lf %lf %lf", atype, &a, &b, &c, &d, &e, &f, &g, &h);
		if( strcmp(atype, atomtype[i])!=0 ) 
			{printf("atom type [%d] is unconsistent! uplimit\n",i); exit(0);}
		uplimit[ipara] = a; ipara++;
		uplimit[ipara] = b; ipara++;
		uplimit[ipara] = c; ipara++;
		uplimit[ipara] = d; ipara++;
		uplimit[ipara] = e; ipara++;
		uplimit[ipara] = f; ipara++;
		uplimit[ipara] = g; ipara++;
		uplimit[ipara] = h; ipara++;
		for(j=1; j<=3; j++){
			fgets(buff, 512, fp);
			sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f, &g, &h);
		uplimit[ipara] = a; ipara++;
		uplimit[ipara] = b; ipara++;
		uplimit[ipara] = c; ipara++;
		uplimit[ipara] = d; ipara++;
		uplimit[ipara] = e; ipara++;
		uplimit[ipara] = f; ipara++;
		uplimit[ipara] = g; ipara++;
		uplimit[ipara] = h; ipara++;
		}
	}
	if(num!=num_ovar) {printf("num_ovar is unconsistent! uplimit\n"); exit(0);}

	//Read the pair-atom parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	//fgets(buff, 512, fp); 
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &a, &b, &c, &d, &e, &f, &g, &h);
		sprintf(name, " %2d %2d",at1,at2); 
		if( strcmp(name, pname[i])!=0 )
			{printf("pair-atom parameters of (%d,%d) are unconsistent! uplimit\n",at1,at2); exit(0);}
		uplimit[ipara] = a; ipara++;
		uplimit[ipara] = b; ipara++;
		uplimit[ipara] = c; ipara++;
		uplimit[ipara] = d; ipara++;
		uplimit[ipara] = e; ipara++;
		uplimit[ipara] = f; ipara++;
		uplimit[ipara] = g; ipara++;
		uplimit[ipara] = h; ipara++;
		fgets(buff, 512, fp);
		sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f, &g, &h);
		uplimit[ipara] = a; ipara++;
		uplimit[ipara] = b; ipara++;
		uplimit[ipara] = c; ipara++;
		uplimit[ipara] = d; ipara++;
		uplimit[ipara] = e; ipara++;
		uplimit[ipara] = f; ipara++;
		uplimit[ipara] = g; ipara++;
		uplimit[ipara] = h; ipara++;
	}	
	if(num!=num_pvar) {printf("num_pvar is unconsistent! uplimit\n"); exit(0);}
		
	//Read the vdw parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %lf %lf %lf %lf %lf %lf", &at1, &at2, &a, &b, &c, &d, &e, &f);
		sprintf(name, " %2d %2d",at1,at2); 
		if( strcmp(name, vdwname[i])!=0 )
			{printf("vdW parameters of (%d,%d) are unconsistent! uplimit\n",at1,at2); exit(0);}
		uplimit[ipara] = a; ipara++;
		uplimit[ipara] = b; ipara++;
		uplimit[ipara] = c; ipara++;
		uplimit[ipara] = d; ipara++;
		uplimit[ipara] = e; ipara++;
		uplimit[ipara] = f; ipara++;
	}
	if(num!=num_vdwvar) {printf("num_vdwvar is unconsistent! uplimit\n"); exit(0);}

	//Read the triple-atom (angle) parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &at3, &a, &b, &c, &d, &e, &f, &g);
		sprintf(name, " %2d %2d %2d",at1,at2,at3); 
		if( strcmp(name, tname[i])!=0 )
			{printf("triple-atom parameters of (%d,%d,%d) are unconsistent! uplimit\n",at1,at2,at3); exit(0);}
		uplimit[ipara] = a; ipara++;
		uplimit[ipara] = b; ipara++;
		uplimit[ipara] = c; ipara++;
		uplimit[ipara] = d; ipara++;
		uplimit[ipara] = e; ipara++;
		uplimit[ipara] = f; ipara++;
		uplimit[ipara] = g; ipara++;
	}	
	if(num!=num_tvar) {printf("num_tvar is unconsistent! uplimit\n"); exit(0);}
		
	//Read the four-atoms (torsion) parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &at3, &at4, &a, &b, &c, &d, &e, &f, &g);
		sprintf(name, " %2d %2d %2d %2d",at1,at2,at3,at4); 
		if( strcmp(name, fname[i])!=0 )
			{printf("4body parameters of (%d,%d,%d,%d) are unconsistent! uplimit\n",at1,at2,at3,at4); exit(0);}
		uplimit[ipara] = a; ipara++;
		uplimit[ipara] = b; ipara++;
		uplimit[ipara] = c; ipara++;
		uplimit[ipara] = d; ipara++;
		uplimit[ipara] = e; ipara++;
		uplimit[ipara] = f; ipara++;
		uplimit[ipara] = g; ipara++;
	}	
	if(num!=num_fvar) {printf("num_fvar is unconsistent! uplimit\n"); exit(0);}
		
	//Read the hydrogen bond parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %lf %lf %lf %lf", &at1, &at2, &at3, &a, &b, &c, &d);
		sprintf(name, " %2d %2d %2d",at1,at2,at3); 
		if( strcmp(name, hname[i])!=0 )
			{printf("hydrogen bond of (%d,%d,%d) are unconsistent! uplimit\n",at1,at2,at3); exit(0);}
		uplimit[ipara] = a; ipara++;
		uplimit[ipara] = b; ipara++;
		uplimit[ipara] = c; ipara++;
		uplimit[ipara] = d; ipara++;
	}	
	if(num!=num_hvar) {printf("num_hvar is unconsistent! uplimit\n"); exit(0);}
	
	fclose(fp);

	if(ipara!=nparas) {printf("Number of total parameters is not consistent to uplimit!\n"); exit(0);}
	return 0;
}


double ReadDownlimit4reaxff(){

	int i, j, ipara, num;
	int at1, at2, at3, at4;
	double a, b, c, d, e, f, g, h;
	char atype[2], name[20];
	char buff[512];
	FILE *fp;	

	ipara=0;

	fp = fopen("downlimit.reax", "r");
	if(!fp){
		printf("There is no downlimit.reax file!!!\n");
		exit(0);
	}
	
	//Read the first line
	fgets(buff, 512, fp);
		
	//Read the global parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%lf %*s", &a);
		downlimit[ipara] = fabs(a); ipara++;
	}
	if(num!=num_gvar) { printf("num_gvar is unconsistent! downlimit\n"); exit(0);}

	//Read the one-atom parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%s %lf %lf %lf %lf %lf %lf %lf %lf", atype, &a, &b, &c, &d, &e, &f, &g, &h);
		if( strcmp(atype, atomtype[i])!=0 ) 
			{printf("atom type [%d] is unconsistent! downlimit\n",i); exit(0);}
		downlimit[ipara]=a; ipara++;
		downlimit[ipara]=b; ipara++;
		downlimit[ipara]=c; ipara++;
		downlimit[ipara]=d; ipara++;
		downlimit[ipara]=e; ipara++;
		downlimit[ipara]=f; ipara++;
		downlimit[ipara]=g; ipara++;
		downlimit[ipara]=h; ipara++;
		for(j=1; j<=3; j++){
			fgets(buff, 512, fp);
			sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f, &g, &h);
		downlimit[ipara]=a; ipara++;
		downlimit[ipara]=b; ipara++;
		downlimit[ipara]=c; ipara++;
		downlimit[ipara]=d; ipara++;
		downlimit[ipara]=e; ipara++;
		downlimit[ipara]=f; ipara++;
		downlimit[ipara]=g; ipara++;
		downlimit[ipara]=h; ipara++;
		}
	}
	if(num!=num_ovar) {printf("num_ovar is unconsistent! downlimit\n"); exit(0);}

	//Read the pair-atom parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &a, &b, &c, &d, &e, &f, &g, &h);
		sprintf(name, " %2d %2d",at1,at2); 
		if( strcmp(name, pname[i])!=0 )
			{printf("pair-atom parameters of (%d,%d) are unconsistent! downlimit\n",at1,at2); exit(0);}
		downlimit[ipara]=a; ipara++;
		downlimit[ipara]=b; ipara++;
		downlimit[ipara]=c; ipara++;
		downlimit[ipara]=d; ipara++;
		downlimit[ipara]=e; ipara++;
		downlimit[ipara]=f; ipara++;
		downlimit[ipara]=g; ipara++;
		downlimit[ipara]=h; ipara++;
		fgets(buff, 512, fp);
		sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f, &g, &h);
		downlimit[ipara]=a; ipara++;
		downlimit[ipara]=b; ipara++;
		downlimit[ipara]=c; ipara++;
		downlimit[ipara]=d; ipara++;
		downlimit[ipara]=e; ipara++;
		downlimit[ipara]=f; ipara++;
		downlimit[ipara]=g; ipara++;
		downlimit[ipara]=h; ipara++;
	}	
	if(num!=num_pvar) {printf("num_pvar is unconsistent! downlimit\n"); exit(0);}
		
	//Read the vdw parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %lf %lf %lf %lf %lf %lf", &at1, &at2, &a, &b, &c, &d, &e, &f);
		sprintf(name, " %2d %2d",at1,at2); 
		if( strcmp(name, vdwname[i])!=0 )
			{printf("vdW parameters of (%d,%d) are unconsistent! downlimit\n",at1,at2); exit(0);}
		downlimit[ipara]=a; ipara++;
		downlimit[ipara]=b; ipara++;
		downlimit[ipara]=c; ipara++;
		downlimit[ipara]=d; ipara++;
		downlimit[ipara]=e; ipara++;
		downlimit[ipara]=f; ipara++;
	}
	if(num!=num_vdwvar) {printf("num_vdwvar is unconsistent! downlimit\n"); exit(0);}

	//Read the triple-atom (angle) parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &at3, &a, &b, &c, &d, &e, &f, &g);
		sprintf(name, " %2d %2d %2d",at1,at2,at3); 
		if( strcmp(name, tname[i])!=0 )
			{printf("triple-atom parameters of (%d,%d,%d) are unconsistent! downlimit\n",at1,at2,at3); exit(0);}
		downlimit[ipara]=a; ipara++;
		downlimit[ipara]=b; ipara++;
		downlimit[ipara]=c; ipara++;
		downlimit[ipara]=d; ipara++;
		downlimit[ipara]=e; ipara++;
		downlimit[ipara]=f; ipara++;
		downlimit[ipara]=g; ipara++;
	}	
	if(num!=num_tvar) {printf("num_tvar is unconsistent! downlimit\n"); exit(0);}
		
	//Read the four-atoms (torsion) parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf", &at1, &at2, &at3, &at4, &a, &b, &c, &d, &e, &f, &g);
		sprintf(name, " %2d %2d %2d %2d",at1,at2,at3,at4); 
		if( strcmp(name, fname[i])!=0 )
			{printf("4body parameters of (%d,%d,%d,%d) are unconsistent! downlimit\n",at1,at2,at3,at4); exit(0);}
		downlimit[ipara]=a; ipara++;
		downlimit[ipara]=b; ipara++;
		downlimit[ipara]=c; ipara++;
		downlimit[ipara]=d; ipara++;
		downlimit[ipara]=e; ipara++;
		downlimit[ipara]=f; ipara++;
		downlimit[ipara]=g; ipara++;
	}	
	if(num!=num_fvar) {printf("num_fvar is unconsistent! downlimit\n"); exit(0);}
		
	//Read the hydrogen bond parameters
	fgets(buff, 512, fp);
	sscanf(buff, "%d %*s", &num);
	for (i=0; i<num; i++){
		fgets(buff, 512, fp);
		sscanf(buff, "%d %d %d %lf %lf %lf %lf", &at1, &at2, &at3, &a, &b, &c, &d);
		sprintf(name, " %2d %2d %2d",at1,at2,at3); 
		if( strcmp(name, hname[i])!=0 )
			{printf("hydrogen bond of (%d,%d,%d) are unconsistent! downlimit\n",at1,at2,at3); exit(0);}
		downlimit[ipara]=a; ipara++;
		downlimit[ipara]=b; ipara++;
		downlimit[ipara]=c; ipara++;
		downlimit[ipara]=d; ipara++;
	}	
	if(num!=num_hvar) {printf("num_hvar is unconsistent! downlimit\n"); exit(0);}
	
	fclose(fp);

	if(ipara!=nparas) {printf("Number of total parameters is not consistent to downlimit!\n"); exit(0);}
	return 0;
}

/* Write data for lammps
Note that atom _ types is fixed to 9, which needs to be modified according to the actual situation.*/
double write_data_md(int obj, int ip){
	int i;
	double box[3];
	double position_x, position_y, position_z;
	FILE *fp;
	box[0] = opt[obj].box[ip][0];
	box[1] = opt[obj].box[ip][1];
	box[2] = opt[obj].box[ip][2];
	
	fp = fopen("data.md", "w");
	fprintf(fp, "# LAMMPS data file\n");
	fprintf(fp, "%d atoms\n",opt[obj].atomicity[ip]);
    fprintf(fp, "%s atom types\n\n","9");
    fprintf(fp,"-10 %f xlo xhi\n",box[0]);
    fprintf(fp,"-10 %f ylo yhi\n",box[1]);
    fprintf(fp,"-10 %f zlo zhi\n",box[2]);
    fprintf(fp,"\nAtoms\n\n");
    for(i=0;i<opt[obj].atomicity[ip];i++){
        fprintf(fp,"%d %d %d %d %f %f %f\n",i+1,1,opt[obj].atom_type[ip][i],0,opt[obj].posx[ip][i],opt[obj].posy[ip][i],opt[obj].posz[ip][i]);
    }
	fclose(fp);
	return 1;
}

/*Wtite input file for lammps*/
double write_in_input(int obj, int ip){
    char linedata[200]={0};
    int i,a;
	FILE *fp;
	FILE *fpw;
	fp = fopen("in.input","r");
	fpw=fopen("in.input2","w");
	while (fgets(linedata,sizeof(linedata)-1,fp))
    {
        if (strncmp(linedata,"pair_coeff",9)==0)
        {
			fprintf(fpw,"%s * * %s ","pair_coeff","para.reax");
            fprintf(fpw,"%s ","C H O N Si S Fe Ni Al");
            fprintf(fpw,"\n");
		}
        else{
            fputs(linedata,fpw);            
    }
	}
    fclose(fp);
    fclose(fpw);
	return(1);
		
}