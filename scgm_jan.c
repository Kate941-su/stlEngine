
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
#define ND (Nx*Ny)
#define NSTEPS 200000
#define NDEVIDE 2000
#define NWRITE (NSTEPS/NDEVIDE)

//prototype function
void gmap_create(int NP, double RX[], double RY[], double L_GX, double L_GY, int N_GX, int N_GY, int NEIGHBOR_ROW[], int NEIGHBOR_COL[], int NEITGHBOR_LEN, int PAIRLIST[][10], double LX, int **G_MAP);

void force(int NP, double RX[], double RY[], double AX[], double AY[], double LX, double LY, int PAIRLIST[][10], double RC2, double* POT, double* POT_IJ, double* MINI);

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
		if (RY[i] > LY - 0.5 || RY[i] < 0.5) {
			VY[i] = -VY[i];
		}
	}
}

void forbid_area(double RX[],double RY[],double RX_B[],double VX[],double VY[],double F_LEFT,double F_RIGHT,double F_TOP,double F_UNDER,int NP){
    int i;
    for(i=0;i<NP;i++){
        //both side
        if((F_TOP+0.5)>RY[i] && (F_UNDER-0.5)<RY[i]){
            if(RX_B[i]>(F_RIGHT+0.5) && (F_RIGHT+0.5)>RX[i]){
                VX[i]=-VX[i];
            }
			else if (RX_B[i]<(F_LEFT-0.5) && RX[i]>(F_LEFT-0.5)){
				VX[i]=-VX[i];
			}
        }
    }
}


void forbid_area2(double RX[],double RY[],double RY_B[],double VX[],double VY[],double F_LEFT,double F_RIGHT,double F_TOP,double F_UNDER,int NP){
    int i;
	for(i=0;i<NP;i++){
        if((F_RIGHT+0.5)>RX[i] && (F_LEFT-0.5)<RX[i]){
			
            if((F_TOP+0.5)<RY_B[i] && (F_TOP+0.5)>RY[i]){
				VY[i]=-VY[i];
				}
        
		    if((F_UNDER-0.5)>RY_B[i] && (F_UNDER-0.5)<RY[i]){
				VY[i]=-VY[i];
				}
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

    FILE *rx_list;
    FILE *ry_list;
    FILE *lxly;
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

	double kin0=0;
	double pot = 0;
	double pot_ij = 0;
	double pot_lis[NWRITE];
	double total_kin=0;
	double total_kin_lis[NWRITE];
	double total_e=0;
	double total_e_lis[NWRITE];


	int i, j, k, l,m;
	int t;
	double ts = 0;
	double te = 0;
	double t_lis[NWRITE];
	double rx_p[ND];
	double ry_p[ND];
	double vx_p[ND];
	double vy_p[ND];
	double ax_p[ND];
	double ay_p[ND];



	//high speed parametor
	const double l_gx = 0.448;
	const double l_gy = l_gx;
	const int n_gx = ceil(lx / l_gx) + 6;//??????????????????
	const int n_gy = ceil(ly / l_gy) + 6;//??????????????????
	const int n_g_all = n_gx * n_gy;
	int pairlist[ND][10];
//	int g_map[230][118];//after change malloc
    int **g_map;
    g_map = (int**)malloc(sizeof(int*)*n_gy);
    for(i=0;i<n_gy;i++){
        g_map[i]=(int*)malloc(sizeof(int)*n_gx);
    }
	int neighbor_list_row[] = { -3,-3,-3,
							   -2,-2,-2,
								-1,-1,-1,
								0,0,0,
								1,1,1,1,
								2,2,2,2,
								3,3,3,3 };

	int neighbor_list_col[] = { -3,-2,-1,
								-3,-2,-1,
								-3,-2,-1,
								-3,-2,-1,
								-3,-2,-1,0,
								-3,-2,-1,0,
								-3,-2,-1,0
	};


	/*
	-3-3-3-3-3-3-3
	-2-2-2-2-2-2-2
	-1-1-1-1-1-1-1
	0 0 0   0 0 0
	1 1 1 1 1 1 1
	2 2 2 2 2 2 2
	3 3 3 3 3 3 3

	*/
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
	for (i = 0;i < ND;i++) {
		for (j = 0;j < 10;j++) {
			pairlist[i][j] = -1;
		}
	}
	
	for (i = 0;i < n_gy;i++) {
		for (j = 0;j < n_gx;j++) {
			g_map[i][j] = -1;
		}
	}
	const int neighbor_len = sizeof(neighbor_list_row) / sizeof(int);
	k = 0;
	for (i = 0;i < Nx;i++) {
		for (j = 0;j < Ny;j++) {
			rx_p[k] = i * dx + dx * 0.5;
			ry_p[k] = j * dy + dy * 0.5;
			vx_p[k] = sqrt(temp) * normal();
			vy_p[k] = sqrt(temp) * normal();
			ax_p[k] = 0.0;
			ay_p[k] = 0.0;
			k = k + 1;
		}
	}


    //12/24 DP?????????????????????????????????
    double f_top=0.5*ly;
    double f_under=0.2*ly;
    double f_left=0.2*lx;
    double f_right=0.8*lx;
    int count_fzone=0;
    int count_buffer[ND];

    for (i=0;i<ND;i++){
        count_buffer[i]=-1;
    }

    for (i=0;i<ND;i++){
        if(f_left-0.5 > rx_p[i] || f_right+0.5 < rx_p[i] ||f_top+0.5 < ry_p[i] || f_under-0.5 > ry_p[i]){
            count_buffer[count_fzone]=i;
            count_fzone+=1;
        }
    }
    double *rx;
    double *vx;
    double *ry;
    double *vy;
    double *ay;
    double *ax;

    double *rx_b;
    double *ry_b;

    int xxx,xx=0;
    rx=(double*)malloc(sizeof(double)*count_fzone);
    ry=(double*)malloc(sizeof(double)*count_fzone);
    vx=(double*)malloc(sizeof(double)*count_fzone);
    vy=(double*)malloc(sizeof(double)*count_fzone);
    ax=(double*)malloc(sizeof(double)*count_fzone);
    ay=(double*)malloc(sizeof(double)*count_fzone);

    rx_b=(double*)malloc(sizeof(double)*count_fzone);
    ry_b=(double*)malloc(sizeof(double)*count_fzone);

    for (i=0;i<count_fzone;i++){
        xxx=count_buffer[i];
        rx[i]=rx_p[xxx];
        ry[i]=ry_p[xxx];
        vx[i]=vx_p[xxx];
        vy[i]=vy_p[xxx];
        ax[i]=ax_p[xxx];
        ay[i]=ay_p[xxx];   
        xx+=1;    
    }

    int N =count_fzone;

    //prepare read file

    char moge[5];
    char vx_name[] = "vx";
    char vy_name[] = "vy";
    double vx_dummy=0.0;
    double vy_dummy=0.0;
    char text[] =".dat";
    char vy_lis[30];
    char vx_lis[30];
    sprintf(moge,"%d",ND);
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

    /////////////////////////////


	for (i = 0;i < N;i++) {
		kin0 = kin0 + (vx[i] * vx[i] + vy[i] * vy[i]);
	}


	kin0 = kin0 * 0.5;

	ts = omp_get_wtime();

	l = 0;
	for (i = 0;i < N;i++)
	{
		for (j = 0;j < 10;j++) {
			if (pairlist[i][j] != -1) {
				l += 1;
			}
		}

	}
	gmap_create(N, rx, ry, l_gx, l_gy, n_gx, n_gy, neighbor_list_row, neighbor_list_col, neighbor_len, pairlist, lx, g_map);
		m=0;
		for (i=3;i<n_gy-3;i++){
			for(j=3;j<n_gx-3;j++){
				if(g_map[i][j]!=-1){
					l+=1;
				}
			}
		}
		printf("particlenum:%d\n",l);
	printf("hello mainloop\n");
	printf("kin0:%lf\n\n",kin0);
tot_file=fopen("./plotdata/totene2.dat","w");
	//-------------start mainroop------------------
	for (t = 1;t <= NSTEPS;t++)
	{
		//initialize
		total_kin = 0.0;
		total_e = 0.0;
//		t_lis[t-1] = t * h;

		//verlet
        for(i=0;i<N;i++){
            rx_b[i]=rx[i];
			ry_b[i]=ry[i];
        }

		for (i = 0;i < N;i++) {
			rx[i] = rx[i] + vx[i] * h + ax[i] * h2;
			ry[i] = ry[i] + vy[i] * h + ay[i] * h2;
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}
	gmap_create(N, rx, ry, l_gx, l_gy, n_gx, n_gy, neighbor_list_row, neighbor_list_col, neighbor_len, pairlist, lx, g_map);
//debug zone
		m=0;
		for (i=3;i<n_gy-3;i++){
			for(j=3;j<n_gx-3;j++){
				if(g_map[i][j]!=-1){
					m+=1;
				}
			}
		}

		if(m!= N){
			printf("particlenum:%d\n",m);
			printf("kasanari\n\n");
			break;
		}

		//cliculate force		
		force(N, rx, ry, ax, ay, lx, ly, pairlist, rc2, &pot, &pot_ij, &mini);

		//second verlet
		for (i = 0;i < N;i++) {
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}

		//boundary condition 
		boundary(N, rx, lx);
		elastic(N, ry, vy, ly);
        forbid_area(rx,ry,rx_b,vx,vy,f_left,f_right,f_top,f_under,N);
        forbid_area2(rx,ry,ry_b,vx,vy,f_left,f_right,f_top,f_under,N);
		//caliclate kin-energy
		for (i = 0;i < N;i++) {
			total_kin += (vx[i] * vx[i] + vy[i] * vy[i]);
		}





		total_kin = total_kin * 0.5;
		total_e = total_kin + pot;

		//write in file
		if (t%10000==0){
			fprintf (tot_file,"%lf\n",total_e);
		}

		if (t % nout == 0)
		{
			printf("Time:%lf,Kinetic Energy:%lf,Potential Energy:%lf,Total Energy:%lf\n", t * h, total_kin, pot, total_e);
			printf("particlenum:%d\n",m);

		}

	}
	//end mainloop
	te = omp_get_wtime();
	fclose(tot_file);
	printf("time cost is %lf seconds\n", te - ts);
	l = 0;

    //free activated array

    for(i=0;i<n_gy;i++){
        free(g_map[i]);
    }
    free(g_map);
    free(rx);
    free(ry);
    free(vx);
    free(vy);
    free(ax);
    free(ay);
    free(rx_b);
    free(ry_b);

//check the gmap
/*
	printf("mini:%lf\n", mini);
	for (i = 0;i < n_gy;i++)
	{
		for (j = 0;j < n_gx;j++) {
			if (g_map[i][j] != -1) {
				l += 1;
			}
		}

	}

	for(i=0;i<N;i++){
		for(j=0;j<10;j++){
			printf("%d  ",pairlist[i][j]);
		}
		printf("\n");
	}
*/

/*
	printf("counter:%d\n", l);
	printf("n_gx:%d\n", n_gx);
	printf("n_gy:%d\n", n_gy);
*/


//file output
/*
kin_file = fopen("./plotdata/kin.dat","w");
pot_file = fopen("./plotdata/pot.dat","w");
tot_file = fopen("./plotdata/tot.dat","w");
enegy_file=fopen("./plotdata/energy.dat","w");
for(i=0;i<NSTEPS;i++){
	if(i%10000==0)
		{
		fprintf(kin_file,"%lf	%lf\n",t_lis[i],total_kin_lis[i]);
		fprintf(pot_file,"%lf	%lf\n",t_lis[i],pot_lis[i]);
		fprintf(tot_file,"%lf	%lf\n",t_lis[i],total_e_lis[i]);
		fprintf(enegy_file,"%lf    %lf    %lf    %lf\n",t_lis[i],total_kin_lis[i],pot_lis[i],total_e_lis[i]);
		}
}
fclose(kin_file);
fclose(pot_file);
fclose(tot_file);
fclose(enegy_file);
fclose(rx_list);
fclose(ry_list);

lxly=fopen("./plotdata/anidat/makefig/lxly.dat","w");
fprintf(lxly,"%lf ",lx);
fclose(lxly);
*/
	return 0;

}






//prototype function after write
void gmap_create(int NP, double RX[], double RY[], double L_GX, double L_GY, int N_GX, int N_GY, int NEIGHBOR_ROW[], int NEIGHBOR_COL[], int NEITGHBOR_LEN, int PAIRLIST[][10], double LX, int **G_MAP) {
	int i, j, k;
	int gx_map=0;
	int gy_map=0;
	int partcle_counter=0;

	for (i = 0;i < N_GY;i++) {
		for (j = 0;j < N_GX;j++) {
			G_MAP[i][j] = -1;
		}
	}

	for (i = 0;i < NP;i++) {
		for (j = 0;j < 10;j++) {
			PAIRLIST[i][j] = -1;
		}
	}

	for (i = 0;i < NP;i++) {
		gx_map = (int)(RX[i] / L_GX) + 3;
		gy_map = (int)(RY[i] / L_GY) + 3;// N_GY-4 - (int)(RY[i] / L_GY)
		G_MAP[gy_map][gx_map] = i;
	}

	for (i = 3;i < N_GY-3;i++) {
		for (j = 0;j < 3;j++) {
			G_MAP[i][j] = G_MAP[i][N_GX - 6 + j];
		}
	}

	for (i = 3;i < N_GY - 3;i++) {
		for (j = 3;j < N_GX - 3;j++) {
			partcle_counter = 0;
			if (G_MAP[i][j] != -1) {
				for (k = 0;k < NEITGHBOR_LEN;k++)
				{
					if (G_MAP[i + NEIGHBOR_ROW[k]][j + NEIGHBOR_COL[k]] != -1)
					{
						PAIRLIST[G_MAP[i][j]][partcle_counter] = G_MAP[i + NEIGHBOR_ROW[k]][j + NEIGHBOR_COL[k]];
						partcle_counter += 1;

					}
				}

				PAIRLIST[G_MAP[i][j]][9] = partcle_counter;
			}
		}
	}



}

void force(int NP, double RX[], double RY[], double AX[], double AY[], double LX, double LY, int PAIRLIST[][10], double RC2, double* POT, double* POT_IJ, double* MINI) {

	int i, j;
	int roop_num=0;
	double rxij=0.0;
	double ryij=0.0;
	double r2=0.0;
	double ir2=0.0, ir6=0.0;
	double fx=0.0, fy=0.0;


	for (i = 0;i < NP;i++) {
		AX[i] = 0.0;
		AY[i] = 0.0;
	}

	*POT = 0.0;
	*POT_IJ = 0.0;

	for (i = 0;i < NP;i++) {
		roop_num = PAIRLIST[i][9];
		for (j = 0;j < roop_num;j++) {
			rxij = RX[i] - RX[PAIRLIST[i][j]];

			if (rxij >= 0.5 * LX) {
				rxij = rxij - LX;
			}

			else if (rxij < -0.5 * LX) {
				rxij = rxij + LX;
			}

			ryij = RY[i] - RY[PAIRLIST[i][j]];

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
			AX[PAIRLIST[i][j]] = AX[PAIRLIST[i][j]] - fx;
			AY[PAIRLIST[i][j]] = AY[PAIRLIST[i][j]] - fy;

			*POT = *POT + *POT_IJ;

		}
	}

}
