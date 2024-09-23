#ifndef _DEFINE_H_
#define _DEFINE_H_

#define min(x,y) (x<y?x:y)

#define MAX_TARGETS	20
#define MAX_DATAPOINTS	50
#define MAX_ATOMICITY	1000
#define MAX_ATOMTYPES	20
#define MAX_PARAMETERS	10000
#define MAX_PARETOS	100

#define EPST		1.0e-3

//#define KB		1.38e-23
//#define KCAL2J	6.9476953e-21
#define KCALM2EV	0.0433634

//Following constants correspond to the ReaxFF parameter file
#define N_GVAR		39
#define N_OVAR		32
#define N_PVAR		16
#define N_VDWVAR	6
#define N_TVAR		7
#define N_FVAR		7
#define N_HVAR		4

#define LMPMPI		"lmp_mpi"
#define OUTPUT		"output"
#define RUNREAXFF	"lmp_mpi > output"

#endif
