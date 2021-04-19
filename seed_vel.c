#include <stdio.h>
#include "MT.h"
#include <string.h>
#include <math.h>
#include <time.h>
#define N (15*15)

double normal(void) {
	double z = sqrt(-2.0 * log(genrand_real3())) * sin(2.0 * M_PI * genrand_real3());
	return z;
}

int main(void){
    FILE *vx_file;
    FILE *vy_file;

    int i;
    int num=N;
    char moge[5];
    char vx_name[] = "vx";
    char vy_name[] = "vy";
    char text[] =".dat";
    char vy_lis[30];
    char vx_lis[30];
	double vx[N];
	double vy[N];
    double temp=2.5;


	//generate random seed from time
	init_genrand((unsigned)time(NULL));

    sprintf(moge,"%d",num);

    sprintf(vx_lis,"%s%s%s",vx_name,moge,text);
    sprintf(vy_lis,"%s%s%s",vy_name,moge,text);


    vx_file=fopen(vx_lis,"w");
    vy_file=fopen(vy_lis,"w");

    for(i=0;i<N;i++){
			vx[i] = sqrt(temp) * normal();
			vy[i] = sqrt(temp) * normal();
            fprintf(vx_file,"%lf\n",vx[i]);
            fprintf(vy_file,"%lf\n",vy[i]);
    }
    fclose(vx_file);
    fclose(vy_file);



    return 0;
}
