#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "deviance.h"
#include "reaxffio.h"
#include "define.h"
#include "common1.h"
#include "common2.h"


double DeviCopy(double *devi_new, double *devi_old){
	int obj;
	for(obj=0; obj<nobjs; obj++)
		devi_old[obj] = devi_new[obj];
	return 0;
}

double DeviCalc(double *deviation,double *origin_devi){
	int obj, ip;
	double var, total_var;
	double **a, b=0.0;
	
	total_var = 0.0;

	/* Calculating the potential energy */
	a = (double **) malloc(sizeof(double *) * nobjs); 
	for (obj=0; obj<nobjs; obj++){
		a[obj] = (double *) malloc(sizeof(double) * opt[obj].npoints_md);
		write_in_input();
		for (ip=0; ip<opt[obj].npoints_md; ip++){
			a[obj][ip] = EnergyCalc_ReaxFF(obj, ip);
		}

		/* Calculating the energy using to compare with reference value */
		/* Firstly, find the datum point (b) and post-analysis the energy (a)
		 * the datum point is different for different energy type, b is
		 * for bond,  the energy of the last point
		 * for angle, the minimum energy, etc... */
		/* For bond, b equals to the energy of longest bond-distance */
		if     (!strncmp(opt[obj].opt_type, "bond", 4))
			b = a[obj][opt[obj].npoints_ref-1];
		/* For EOS, b = 0, and a should be divided by the atomicity */
		else if( (!strncmp(opt[obj].opt_type, "EOS", 3))||(!strncmp(opt[obj].opt_type, "eos", 3)) ){
			b = 0;
			for(ip=0; ip<opt[obj].npoints_ref; ip++)
				a[obj][ip] /= opt[obj].atomicity[ip];	
		}
		/* For angle and torsion, b should be the minimum energy */
		else if(!strncmp(opt[obj].opt_type, "angle", 5)){
			b = a[obj][0];	
			for(ip=0; ip<opt[obj].npoints_ref; ip++)
				if(a[obj][ip]<b) 
					b = a[obj][ip];
		}
		else if(!strncmp(opt[obj].opt_type, "torsion", 7)){
			b = a[obj][0];
			for(ip=0; ip<opt[obj].npoints_ref; ip++)
				if(a[obj][ip]<b)
					b = a[obj][ip];
		}

		/* Secondly, minus the datum point (b) to get a value whihc is comparable to the reference */
		for(ip=0; ip<opt[obj].npoints_ref; ip++)
			a[obj][ip] -= b;

		/* For other,  */
		//if(!strncmp(opt[obj].opt_type, "other", 5))
		//	SelfDefine_Calculation(a, obj);

		/* For charge */
		/*else if(!strcmp(opt[obj].opt_type, "charge")){
		}*/

		/* Secondly, minus the datum point (b) */
		/* and calculate the deviance to ref.  */
		deviation[obj] = 0.0;
		origin_devi[obj] = 0.0;
		for(ip=0; ip<opt[obj].npoints_ref; ip++){
			var = a[obj][ip] - opt[obj].ref_axisy[ip];
			deviation[obj] += opt[obj].eta[ip]*var*var;
			origin_devi[obj] += var*var;
		}
		total_var += lamda[obj]*deviation[obj];
	}
	return total_var;
}

double EnergyCalc_ReaxFF(int obj, int ip){

	double a;
	write_reaxff_input_files(obj, ip);
	system("lmp_mpi <in.input2");
	sleep(0.1);
	system("rm -f *lammps.e*");
	a = analysis_for_reaxff_output(outname);
	a *= KCALM2EV;
	/* Post-analysis for different energy type */
	/* bond, angle, EOS, torsion, charge, etc..*/

	return a;
}


