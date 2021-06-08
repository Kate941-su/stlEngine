//5/27

//This file is able to use only remote
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "MT.h"
#include <time.h>

#define Nx 2
#define Ny 1
#define N (Nx*Ny)
#define NSTEPS 100000
#define NDEVIDE 1000
#define NWRITE (NSTEPS/NDEVIDE)

//prototype function
void gmap_create(int NP, double RX[], double RY[], double L_GX, double L_GY, int N_GX, int N_GY, int NEIGHBOR_ROW[], int NEIGHBOR_COL[], int NEITGHBOR_LEN, int PAIRLIST[][10], double LX, int **G_MAP);

void force(int NP, double RX[], double RY[], double AX[], double AY[], double LX, double LY, int PAIRLIST[][10], int PAIRLIST2[][10] ,double RC2, double* POT, double* POT_IJ, double* MINI,FILE *PAIR);

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

void elastic(int NP, double RY[], double VY[], double LY,double *QIN,double *QOUT) {
	int i;
	double dq,vl;
	for (i = 0;i < NP;i++) {
		if (RY[i] > LY - 0.5 | RY[i] < 0.5) {
			vl=VY[i];
			VY[i] = -1.0*VY[i];
			dq=0.5*(VY[i]*VY[i]-vl*vl);
			if(dq>0){
				*QIN+=dq;
			} 
			else{
				*QOUT+=dq;
			}
		}
	}
}

bool heatwall(int NP, double RY[], double VY[], double LY,double TEMP_L,double TEMP_H,double *Q_IN,double *Q_OUT){
	int i;
	double dq=0.0;
	double vy_l=0.0;
	for (i = 0;i < NP;i++) {
		if (RY[i] > LY - 0.5) {
            vy_l=VY[i];
			VY[i] = -canon_b(TEMP_H);
			dq=0.5*(VY[i]*VY[i]-vy_l*vy_l);

			if(dq>=0.0){
				*Q_IN+=dq;

			}else{
				*Q_OUT+=dq;
			}
		}

		else if (RY[i] < 0.5) {
            vy_l=VY[i];
			VY[i] = canon_b(TEMP_L);
			dq=0.5*(VY[i]*VY[i]-vy_l*vy_l);
			if(dq>=0.0){
				*Q_IN+=dq;
			}else{
				*Q_OUT+=dq;
			}
		}
	}
    return dq == 0.0 ? true : false;
}

//粒子位置をファイルに書き込む
void writeParticlePos(FILE* rx_list, FILE* ry_list, double* rx, double* ry) {
    for (int i = 0; i < N; i++) {
        fprintf(rx_list , "%lf  ", rx[i]);
        fprintf(ry_list , "%lf  ", ry[i]);
    }
    fprintf(rx_list, "\n");
    fprintf(ry_list, "\n");
}

int main(void) {
	//FILEopen
	FILE *kin_file;
	FILE *pot_file;
	FILE *tot_file;
	FILE *enegy_file;
    FILE *vx_read;
    FILE *vy_read;
	FILE *pair;
	//for animation
	FILE *rx_list;
	FILE *ry_list;
	FILE *ay_list;
	FILE *ax_list;
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
	const double h = 3e-4;
	const double h2 = 0.5 * h * h;
	const int nout = NSTEPS / 10;
	const int t_end = NSTEPS * h;
	int t_write=0;



	double kin0=0;
	double diff=0;
	double diff_tot=0;
	double total_e_b=0;
	double pot = 0;
	double pot_ij = 0;
	double total_kin=0;
	double total_e=0;



	int i, j, k, l,m;
	int t;
	double ts = 0;
	double te = 0;

	double rx[N];
	double ry[N];
	double vx[N];
	double vy[N];
	double ax[N];
	double ay[N];


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


	//high speed parametor
	//ぴったりにする
	int bo=(lx/0.4);
	double l_gx_d=lx/bo;

	const double l_gx = l_gx_d;
	const double l_gy = l_gx;
	const int n_gx = ceil(lx / l_gx) + 6;//袖領域のため
	const int n_gy = ceil(ly / l_gy) + 40;//袖領域のため
	const int n_g_all = n_gx * n_gy;


	int pairlist[N][10];
	int pairlist2[N][10];

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

    double q_in=0;
    double q_out=0;
	
	for (i = 0;i < N;i++) {
		for (j = 0;j < 10;j++) {
			pairlist[i][j] = -1;
		}
	}
	
	for (i = 0;i < N;i++) {
		for (j = 0;j < 10;j++) {
			pairlist2[i][j] = -1;
		}
	}


	for (i = 0;i < n_gy;i++) {
		for (j = 0;j < n_gx;j++) {
			g_map[i][j] = -1;
		}
	}

	const int neighbor_len = sizeof(neighbor_list_row) / sizeof(int);


///read vx,vy,file

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
//			vx[k] = sqrt(temp) * normal();
//			vy[k] = sqrt(temp) * normal();
			ax[k] = 0.0;
			ay[k] = 0.0;
			k = k + 1;
		}
	}


	for (i = 0;i < N;i++) {
		kin0 = kin0 + (vx[i] * vx[i] + vy[i] * vy[i]);
	}


	kin0 = kin0 * 0.5;
	total_e=kin0;
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
    printf("kin0:%lf\n",kin0);
    tot_file=fopen("./plotdata/totgmap2.dat","w");
    FILE* box = fopen("./plotdata/box.dat","w");
    rx_list=fopen("./plotdata/rx_list.dat","w");
    ry_list=fopen("./plotdata/ry_list.dat","w");
    ax_list=fopen("./plotdata/ax_list.dat","w");
    ay_list=fopen("./plotdata/ay_list.dat","w");
    pair=fopen("./plotdata/pairlist_kai.dat","w");

    //シリンダの大きさをファイルに書き込んでおく
    fprintf(box,"%lf    %lf",lx , ly);
    fclose(box);
 //   bool isForce = false;
//    bool isHit = false;
	//-------------start mainroop------------------

	for (t = 1;t <= NSTEPS;t++)
	{
//        if (pot > 0) {
//           isForce = true; 
//        }
        total_e_b=total_e;
		//initialize
		total_kin = 0.0;
		total_e = 0.0;

		if ((t-1)%NDEVIDE==0){
			t_write+=1;
		}
//        isHit = heatwall(N, ry,vy,ly,temp_l,temp_h,&q_in,&q_out);//熱壁条件
        elastic(N ,ry, vy, ly, &q_in, &q_out);
		//verlet
//        #pragma omp parallel for
		for (i = 0;i < N;i++) {
			rx[i] = rx[i] + vx[i] * h + ax[i] * h2;
			ry[i] = ry[i] + vy[i] * h + ay[i] * h2;
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}
//make gmap and pairlist
        boundary(N, rx, lx);
		gmap_create(N, rx, ry, l_gx, l_gy, n_gx, n_gy, neighbor_list_row, neighbor_list_col, neighbor_len, pairlist, lx, g_map);
		//cliculate force		
		force(N, rx, ry, ax, ay, lx, ly, pairlist, pairlist2 , rc2, &pot, &pot_ij, &mini,pair);

		//second verlet
		for (i = 0;i < N;i++) {
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}
		//caliclate kin-energy
		for (i = 0;i < N;i++) {
			total_kin += (vx[i] * vx[i] + vy[i] * vy[i]);
		}
		total_kin = total_kin * 0.5;
		total_e = total_kin + pot;
		diff=total_e-total_e_b;
		diff_tot+=diff;
	//write in file
    if (t % NDEVIDE == 0) {
        fprintf(tot_file , "%lf    %lf    %lf\n" , t * h , kin0, total_e - q_in - q_out);
        writeParticlePos(rx_list, ry_list, rx, ry);
    }

    double stat = diff_tot - q_in - q_out;//デバッグ 
//    if (isHit && isForce) {
//        printf("%lf   %lf\n",q_in + q_out + kin0,total_e);
//    }
	if (t % nout == 0)
		{
        printf("Time:%lf,Kinetic Energy:%lf,Potential Energy:%lf,Total Energy:%lf\n", t * h, total_kin, pot, total_e);
        printf("e stat: %lf\n",(diff_tot-q_in-q_out));
		}
//    isHit = false;
//    isForce = false;
	}

//end mainloop
//close particle animation
	fclose(rx_list);
	fclose(ry_list);
	fclose(ax_list);
	fclose(ay_list);
	fclose(tot_file);
	fclose(pair);
	te = omp_get_wtime();
//	printf("time cost is %lf seconds\n", te - ts);
	l = 0;
    //free activated array
    for(i=0;i<n_gy;i++){
        free(g_map[i]);
    }
    free(g_map);
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
		gy_map =N_GY- 4 - (int)(RY[i] / L_GY);// (int)(RY[i] / L_GY) + 3
		G_MAP[gy_map][gx_map] = i;
	}

	for (i = 0;i < N_GY;i++) {
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

void force(int NP, double RX[], double RY[], double AX[], double AY[], double LX, double LY, int PAIRLIST[][10], int PAIRLIST2[][10] ,double RC2, double* POT, double* POT_IJ, double* MINI,FILE *PAIR) {

	int i, j,k=0;
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
		k=0;
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

				PAIRLIST2[i][k]=PAIRLIST[i][j];
			}

			AX[i] = AX[i] + fx;
			AY[i] = AY[i] + fy;
			AX[PAIRLIST[i][j]] = AX[PAIRLIST[i][j]] - fx;
			AY[PAIRLIST[i][j]] = AY[PAIRLIST[i][j]] - fy;

			*POT = *POT + *POT_IJ;

		}
	}

}
