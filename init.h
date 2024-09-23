#ifndef _INIT_H_
#define _INIT_H_

#include "deviance.h"

int InitSys();
int InitiationPartOne();
int InitiationPartTwo();
int InitiationPSO();
int ReadFiles();
int ReadInputInfo();
int LamdaAssign();
int InitCalc();

double find_atom_weight(char type[2]);
char   find_atom_type(char an[2]);
#endif

