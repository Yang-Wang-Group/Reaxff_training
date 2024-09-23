#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <math.h>

#include "define.h"
#include "common1.h"
#include "vardefi.h"
#include "common2.h"

#include "init.h"
#include "deviance.h"
#include "mcmosa.h"
#include "output.h"
#include "reaxffio.h"
#include "pso.h"

int main(){
	//MC-MOSA method
	
	int flag_accept;
	int count;
	double possibility;

	InitSys();

	/* input information */
	temp = InitTemp;
	while (temp>FinalTemp){
		flag_accept = 0;
		count = 0;
		while (!flag_accept){
			RandomMove();
			DeviCalc(devi_new,origin_deviation);
			possibility = McmosaCompare(devi_new, devi_old);
			flag_accept = McmosaUpdata(possibility, devi_new);
			count++;
			/* If the current pareto result is diffuclt to be accepted, then make *
			 * the current pareto to be the best pareto and optimize it again     */
			if(count>=interval2){
				FindBestPareto(1);
				count = 0;
			}
		}
		WriteOutput(1);
		CoolTemp();

		nloop+= 1;
	}

	WriteOutputParetos();
	return 1;
}


