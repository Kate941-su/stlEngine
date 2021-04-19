//g_mapの配列をmallocすることでもっといいコードがかけるかもしれないがc++で実装するのもあり
//とりあえず静的な二次元配列を確保しておく


//This file is able to use only remote
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "MT.h"
#include <time.h>

#define Nx 10
#define Ny 20
#define N (Nx*Ny)
#define NSTEPS 10000


//prototype function

void force(int NP, double RX[], double RY[], double AX[], double AY[], double LX, double LY,double RC2, double* POT, double* POT_IJ, double* MINI);

void bk(int NP,double RX[],double RY[],double RC, double LX,int PAIRLIST[][10]);
//normal distribution function
//genrand_real3 is imported form MT.h
double normal(void) {
	double z = sqrt(-2.0 * log(genrand_real3())) * sin(2.0 * M_PI * genrand_real3());
	return z;
}

double canon_a(double TEMP) {
	double res;
	res = sqrt(TEMP) * normal();
	return res;
}

double canon_b(double TEMP) {
	double res;
	res = sqrt(-2.0 * TEMP * log(genrand_real3()));
	return res;
}
/*bondary condition*/
void boundary(int NP, double RX[], double LX) {
	int i;
	for (i = 0;i < NP;i++) {
		if (RX[i] > LX) {
			RX[i] = RX[i] - LX;
		}
		else if (RX[i] <= 0.0) {
			RX[i] = RX[i] + LX;
		}
	}
}

void elastic(int NP, double RY[], double VY[], double LY) {
	int i;
	for (i = 0;i < NP;i++) {
		if (RY[i] > LY - 0.5 | RY[i] < 0.5) {
			VY[i] = -VY[i];
		}
	}
}




int main(void) {
	//FILEopen
	FILE *kin_file;
	FILE *pot_file;
	FILE *tot_file;
	FILE *enegy_file;
    FILE *vx_read;
    FILE *vy_read;
	//generate random seed from time
	init_genrand((unsigned)time(NULL));

	//static number
	const double rho = 0.2;
	double dx = sqrt(1.0 / rho);
	const double dy = dx;
	const double lx = dx * Nx;
	const double ly = dy * Ny;
	const double temp = 2.5;
	const double temp_l = 10.0;
	const double temp_h = 1.0;
	const double temp_p = (temp_h + temp_l) * 0.5;
	const double h = 1e-3;
	const double h2 = 0.5 * h * h;
	const int nout = NSTEPS / 10;
	const int t_end = NSTEPS * h;

	double q_in=0;
	double q_out = 0;
	double q_in_sum = 0;
	double q_out_sum = 0;
	double q_in_lis[NSTEPS];
	double q_out_lis[NSTEPS];


	double kin0=0;
	double pot = 0;
	double pot_ij = 0;
	double pot_lis[NSTEPS];
	double total_kin=0;
	double total_kin_lis[NSTEPS];
	double total_e=0;
	double total_e_lis[NSTEPS];


	int i, j, k, l,m;
	int t;
	double ts = 0;
	double te = 0;
	double t_lis[NSTEPS];
	double rx[N];
	double ry[N];
	double vx[N];
	double vy[N];
	double ax[N];
	double ay[N];

    int pairlist[N][10];
    for (i=0;i<N;i++){
        for (j=0;j<10;j++){
            pairlist[i][j]=-1;
        }
    }


	//idealy thermal efficiency
	const double eta = 1.0 - (temp_h / temp_l);
	const double eta_ca = 1.0 - (sqrt(temp_h / temp_l));

	//force parametor
	const double rc = pow(2.0, (double)1 / 6);
	double rc2 = rc * rc;
	double mini = 1000.0;

	//debug parameter
	double FX=0.0;
	double FY=0.0;



    //prepare read file
    int num=N;
    char moge[5];
    char vx_name[] = "vx";
    char vy_name[] = "vy";
    double vx_dummy=0.0;
    double vy_dummy=0.0;
    char text[] =".dat";
    char vy_lis[30];
    char vx_lis[30];
    sprintf(moge,"%d",num);
    sprintf(vx_lis,"%s%s%s",vx_name,moge,text);
    sprintf(vy_lis,"%s%s%s",vy_name,moge,text);

    vx_read=fopen(vx_lis,"r");
    vy_read=fopen(vy_lis,"r");

    for (i=0;i<N;i++){
        fscanf(vx_read,"%lf",&vx_dummy);
        fscanf(vy_read,"%lf",&vy_dummy);
        vx[i]=vx_dummy;
        vy[i]=vy_dummy;
    }

    fclose(vx_read);
    fclose(vy_read);

	k = 0;
	for (i = 0;i < Nx;i++) {
		for (j = 0;j < Ny;j++) {
			rx[k] = i * dx + dx * 0.5;
			ry[k] = j * dy + dy * 0.5;
	//		vx[k] = sqrt(temp) * normal();
	//		vy[k] = sqrt(temp) * normal();
			ax[k] = 0.0;
			ay[k] = 0.0;
			k = k + 1;
		}
	}


	for (i = 0;i < N;i++) {
		kin0 = kin0 + (vx[i] * vx[i] + vy[i] * vy[i]);
	}


	kin0 = kin0 * 0.5;

	ts = omp_get_wtime();

//	printf("kin0:%lf\n", kin0);
//	printf("counter:%d\n", l);
	printf("hello mainloop\n");
	printf("kin0:%lf\n",kin0);
	//-------------start mainroop------------------

	for (t = 1;t <= NSTEPS;t++)
	{

		//initialize
		total_kin = 0.0;
		total_e = 0.0;
		t_lis[t-1] = t * h;

		//verlet
		for (i = 0;i < N;i++) {
			rx[i] = rx[i] + vx[i] * h + ax[i] * h2;
			ry[i] = ry[i] + vy[i] * h + ay[i] * h2;
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}


		//cliculate force		
		force(N, rx, ry, ax, ay, lx, ly, rc2, &pot, &pot_ij, &mini);

		//second verlet
		for (i = 0;i < N;i++) {
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}

		//boundary condition 
		boundary(N, rx, lx);

		elastic(N, ry, vy, ly);

		//caliclate kin-energy
		for (i = 0;i < N;i++) {
			total_kin += (vx[i] * vx[i] + vy[i] * vy[i]);
		}

		total_kin = total_kin * 0.5;

		total_kin_lis[t-1] = total_kin;
		pot_lis[t-1] = pot;
		total_e = total_kin + pot;
		total_e_lis[t-1] = total_e;

		if (t % nout == 0)
		{
			printf("Time:%lf,Kinetic Energy:%lf,Potential Energy:%lf,Total Energy:%lf\n", t * h, total_kin, pot, total_e);

		}


		if(fabs(total_kin-kin0)>10){
					break;

		}


	}
	//end mainloop
	te = omp_get_wtime();

	printf("time cost is %lf seconds\n", te - ts);
	l = 0;




//file write
kin_file = fopen("./plotdata/kin.dat","w");
pot_file = fopen("./plotdata/pot.dat","w");
tot_file = fopen("./plotdata/tot.dat","w");
enegy_file=fopen("./plotdata/energy.dat","w");
for(i=0;i<NSTEPS;i++){
	fprintf(kin_file,"%lf	%lf\n",t_lis[i],total_kin_lis[i]);
	fprintf(pot_file,"%lf	%lf\n",t_lis[i],pot_lis[i]);
	fprintf(tot_file,"%lf	%lf\n",t_lis[i],total_e_lis[i]);
	fprintf(enegy_file,"%lf    %lf    %lf    %lf\n",t_lis[i],total_kin_lis[i],pot_lis[i],total_e_lis[i]);
}
fclose(kin_file);
fclose(pot_file);
fclose(tot_file);
fclose(enegy_file);

	return 0;

}






void force(int NP, double RX[], double RY[], double AX[], double AY[], double LX, double LY, double RC2, double* POT, double* POT_IJ, double* MINI) {

	int i, j;
	double rxij;
	double ryij;
	double r2;
	double ir2, ir6;
	double fx, fy;


	for (i = 0;i < NP;i++) {
		AX[i] = 0.0;
		AY[i] = 0.0;
	}

	*POT = 0.0;
	*POT_IJ = 0.0;

	for (i = 0;i < NP-1;i++) {
		for (j = i+1;j < NP;j++) {
			rxij = RX[i] - RX[j];

			if (rxij >= 0.5 * LX) {
				rxij = rxij - LX;
			}

			else if (rxij < -0.5 * LX) {
				rxij = rxij + LX;
			}

			ryij = RY[i] - RY[j];

			r2 = rxij * rxij + ryij * ryij;

			ir2 = 1.0 / r2;
			ir6 = ir2 * ir2 * ir2;

			if (r2 >= RC2) {
				fx = 0.0;
				fy = 0.0;
				*POT_IJ = 0.0;
			}else{
				fx = 24.0 * ir6 * (2.0 * ir6 - 1.0) * ir2 * rxij;
				fy = 24.0 * ir6 * (2.0 * ir6 - 1.0) * ir2 * ryij;
				*POT_IJ = (4.0 * (ir6 * ir6 - ir6) + 1.0);

				if (*MINI > sqrt(r2)) {
					*MINI = sqrt(r2);
				}
			}

			AX[i] = AX[i] + fx;
			AY[i] = AY[i] + fy;
            AX[j] = AX[j] - fx;
			AY[j] = AY[j] - fy;
			*POT = *POT + *POT_IJ;

		}
	}

}

void bk(int NP,double RX[],double RY[],double RC, double LX, int PAIRLIST[][10]){
    int i,j,k=0;
    double rxij,ryij,r;
    for (i=0;i<NP-1;i++){
        k=0;
        for(j=i+1;j<NP;j++){

            rxij=RX[i]-RX[j];
			if (rxij >= 0.5 * LX) {
				rxij = rxij - LX;
			}

			else if (rxij < -0.5 * LX) {
				rxij = rxij + LX;
			}
            ryij=RY[i]-RY[j];

            r=sqrt(rxij*rxij+ryij*ryij);
            if(r<RC){
                PAIRLIST[i][k]=j;
                k+=1;
            }
        }

    }
}