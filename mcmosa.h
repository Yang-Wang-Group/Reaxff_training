#ifndef _MCMOSA_H_
#define _MCMOSA_H_
//double RandomMove(); 
double McmosaCompare(double *devi_new, double *devi_old);
int    McmosaUpdata(double p, double *devi);
int    UpdateParetos(double *devi);
int    DeleteBadPareto();
int    FindBestPareto(int option);
int    CoolTemp();
#endif
