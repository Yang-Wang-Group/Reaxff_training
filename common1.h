#ifndef _COMMON1_H_
#define _COMMON1_H_

#include "define.h"

typedef struct tagOptTargets{
	char opt_type[30];	//bond, angle, EOS, torsion, charge, etc...
	char opt_detail[40];    //H-H, C-H, diamond, etc...
	double ref_axisx[MAX_DATAPOINTS];	//ref_axisx[points]: points means how many points are calculated in this energy curve
	double ref_axisy[MAX_DATAPOINTS];	// 					 and each point is a molecular structure
	double min_energy_axisy;
	double min_energy_axisx;
	double eta[MAX_DATAPOINTS];
	int npoints_ref;
	int npoints_md;
	int atomicity[MAX_DATAPOINTS];		//atomicity[points]: atomicity in each point
	int n_atom_species[MAX_DATAPOINTS];
	int atom_type[MAX_DATAPOINTS][MAX_ATOMICITY];//atom_type[points][id]: 1C, 2H, 3O, 4xxxx
	double posx[MAX_DATAPOINTS][MAX_ATOMICITY];		//posx[points][id];  fraction coordination, id meands the No. of atoms
	double posy[MAX_DATAPOINTS][MAX_ATOMICITY];
	double posz[MAX_DATAPOINTS][MAX_ATOMICITY];
	double box[MAX_DATAPOINTS][3];		//box of the cell
} OptTargets;

typedef struct tagPareto{
	double para[MAX_PARAMETERS];
	double devi[MAX_TARGETS];
	double total_devi;
	int   deleted_remark;
} Pareto;

#endif
