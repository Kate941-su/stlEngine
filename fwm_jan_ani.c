#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "MT.h"
#include <time.h>

#define Nx 15
#define Ny 30
#define PARTNUM 2
#define NB (Nx*Ny)
#define NSTEPS 2000
#define NDEVIDE 100
#define NWRITE (NSTEPS/NDEVIDE)
#define TEMPL 4.0
#define PROBABIRITY 1.1
#define SIZE_DP 0.0

//prototype function
void gmap_create(int NP, double RX[], double RY[], double L_GX, double L_GY, int N_GX, int N_GY, int NEIGHBOR_ROW[], int NEIGHBOR_COL[], int NEITGHBOR_LEN, int PAIRLIST[][10], double LX, int **G_MAP);

void force(int NP, double RX[], double RY[], double AX[], double AY[], double LX, double LY, int PAIRLIST[][10], double RC2, double* POT, double* POT_IJ, double* MINI);

void heatwall(double H,int NP,double RY[],double RY_B[],double VY[],double *Q_IN,double *Q_OUT ,double PPY,double PPV,double TEMP_L,double TEMP_H,double *FPP,double LY,double *h_ss,double *d_w);

void piston_move_u(int NP,double RY[],double RY_B[],double VY[],double VY_B[],double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,double *H1_D,double Momentum_u[2],double Kin_u[2],int *k,int *j,double MDP,double *delta_mom,double *delta_kin,double F_RIGHT,double F_LEFT,double RX[]);

void piston_move_d(int NP,double RY[],double RY_B[],double VY[],double VY_B[],double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,double *H1_D,double Momentum_d[2],double Kin_d[2],int *kk,int *jj,double MDP,double *delta_mom,double *delta_kin,double F_RIGHT,double F_LEFT,double RX[]);

void velocity(int NP,double RY[],double VX[],double VY[],double DPY,double *TEMP_D,double *TEMP_U,int *ND,int*NU ,int N_U_list[],double *TEMP);

void partition(double RY[],double VX[] , double VY[] ,double PPY,double DPY,double TEMP_PART_LIST_DOWN[NWRITE][PARTNUM/2],double TEMP_PART_LIST_UP[NWRITE][PARTNUM/2],double PRESS_PART_LIST_DOWN[NWRITE][PARTNUM/2],double PRESS_PART_LIST_UP[NWRITE][PARTNUM/2],int T_WRITE,int NP);


//return max
double judge_max(double array[],int len_array,int *u_num){
	int i;
	double max=0.0;
	for(i=0;i<len_array;i++){
		if (array[i]>max){
			max=array[i];
		}
	}
	return max;
}


double judge_mini(double array[],int len_array,int *l_num){
	int i;
	double mini=1000.0;
	for (i=0;i<len_array;i++){
		if(array[i]<mini){
			mini=array[i];
			*l_num=i;
		}
	}
	return mini;

}

//sumnation array index

int int_array_sum(int NP,int array[NP]){
//    int len;
    int i,N;
    int sum=0;
    for (i=0;i<NP;i++){
        sum+=array[i];
    }
    return sum;
}

double double_array_sum(int NP,double array[NP]){
    int i,N;
    double sum=0;
    for (i=0;i<NP;i++){
        sum+=array[i];
    }
    return sum;
}

int bigger(int a,int b,int *UP_DOWN_FLAG){
	*UP_DOWN_FLAG=0;
	if(a>b){
		*UP_DOWN_FLAG=1;
		printf("up->down");
		return a;
	}else{
		printf("down->up");
		return b;
	}
};


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

void forbid_area(double RX[],double RY[],double RX_B[],double VX[],double VY[],double F_LEFT,double F_RIGHT,double F_TOP,double F_UNDER,int NP){
    int i;
	double RX_D;
    for(i=0;i<NP;i++){
        //both side
		RX_D=fabs(RX[i]-RX_B[i]);
        if((F_TOP+0.5)>RY[i] && (F_UNDER-0.5)<RY[i] && RX_D<1.0){
            if(RX_B[i]>(F_RIGHT+0.5) && (F_RIGHT+0.5>RX[i])){
                VX[i]=-VX[i];
            }
			else if (RX_B[i]<(F_LEFT-0.5) && RX[i]>(F_LEFT-0.5) && RX_D<1.0){
				VX[i]=-VX[i];
			}
        }
    }
}


void forbid_area2(double RX[],double RY[],double RY_B[],double VX[],double VY[],double F_LEFT,double F_RIGHT,double F_TOP,double F_UNDER,int NP){
    int i;
	for(i=0;i<NP;i++){
        if(F_RIGHT+0.5>RX[i] && F_LEFT-0.5<RX[i]){
			
            if(F_TOP+0.5<RY_B[i] && F_TOP+0.5>RY[i]){
				VY[i]=-VY[i];
				}
        
		    if(F_UNDER-0.5>RY_B[i] && F_UNDER-0.5<RY[i]){
				VY[i]=-VY[i];
				}
        }
	}
}

void fix(double RX[],double RY[],double VX[],double VY[],double F_LEFT,double F_RIGHT,double F_TOP,double F_UNDER,int NP){
	int i;
	double temp_x,temp_y;
	for (i=0;i<NP;i++){
		if (F_RIGHT>RX[i] && F_LEFT<RX[i] && F_TOP>RY[i] && F_UNDER<RY[i]){
			temp_x=VX[i];
			temp_y=VY[i];
			VX[i]=-1*fabs(temp_x);
			VY[i]=-1*fabs(temp_y);
		}
	}

}

int main(void) {
	//FILEopen
	FILE *kin_file;
	FILE *pot_file;
	FILE *tot_file;
	FILE *enegy_file;
	FILE *dis_p;
	FILE *pow_p;
	FILE *omega_file;
	FILE *rx_file;
	FILE *rx_b_file;
	FILE *ry_b_file;
	FILE *ry_file;
	FILE *dp_b;
	FILE *dp;
	FILE *py;
	FILE *lxly;
    FILE *vx_read;
    FILE *vy_read;
	FILE *press_u_file;
	FILE *press_d_file;
	//for animation
	FILE *rx_list;
	FILE *ry_list;
	FILE *theta_file;
    FILE *efile;
    FILE *dpx_file;
    FILE *pwx_file;
    FILE *dwx_file;
    FILE *pwy_file;
    FILE *dwy_file;
    FILE *rf_file;


	FILE *press_file;
	FILE *temp_file;
	FILE *part_press_file;
	FILE *part_temp_file;
	FILE *all_condition_file;

	//flexible dir name
	char name_kin[50];
	char name_pot[50];
	char name_tot[50];
	char name_ene[50];
	char name_dis[50];
	char name_pow[50];
	char name_ome[50];
	char name_the[50];
	char name_e[50];
	char name_press_d[50];
	char name_press_u[50];
	char name_press[50];
	char name_temp[50];
	char name_press_part[50];
	char name_temp_part[50];
	char name_all_condition[50];
/*
	int A=10*PROBABIRITY;
	switch (A)
	{
	case 1:
		char plotdata[20]="./plotdata/pro1temp";		
		break;
	case 2:
		char plotdata[20]="./plotdata/pro2temp";		
		break;
	case 3:
		char plotdata[20]="./plotdata/pro3temp";		
		break;
	case 4:
		char plotdata[20]="./plotdata/pro4temp";		
		break;
	default:
		char plotdata[20]="./plotdata/pro5temp";		
		break;
	}
*/

	char plotdata[20]="./plotdata/ani5temp";
	char c_kin[20]="/kin.dat";
	char c_pot[20]="/pot.dat";
	char c_tot[20]="/tot.dat";
	char c_ene[20]="/energy.dat";
	char c_dis[20]="/dis_p.dat";
	char c_pow[20]="/pow_p.dat";
	char c_ome[20]="/omega.dat";
	char c_the[20]="/theta.dat";
	char c_e[20]="/e_lis.dat";
	char c_press_d[20]="/press_d.dat";
	char c_press_u[20]="/press_u.dat";
	char c_press[20]="/press.dat";
	char c_temp[20]="/temp.dat";
	char c_press_part[20]="/press_part.dat";
	char c_temp_part[20]="/temp_part.dat";
	char c_all_condition[20]="/all_condition.dat";
	//generate random seed from time
	init_genrand((unsigned)time(NULL));

	//static number
	const double rho = 0.2;
	double dx = sqrt(1.0 / rho);
	const double dy = dx;
	const double lx = dx * Nx;
	const double ly = dy * Ny;
	const double temp = 2.5;
	const double temp_l = TEMPL;
	const double temp_h = 1.0;
	const double temp_p = (temp_h + temp_l) * 0.5;
	const double h = 1e-3;
	const double h_rev=1/h;
	const double h2 = 0.5 * h * h;
	const int nout = NSTEPS / 10;
	const int t_end = NSTEPS * h;
	
	double kin0=0;
	double pot = 0;
	double pot_ij = 0;
	double total_kin=0;
	double total_e=0;
	const int pow_num=4;
//	check momentum and energy
	double momentum_u[2]={0.0,0.0};
	double kin_u[2]={0.0,0.0};
	double momentum_d[2]={0.0,0.0};
	double kin_d[2]={0.0,0.0};

	int up_hit=0;
	int down_hit=0;
    int up_through=0;
    int down_through=0;
	double delta_mom_u=0.0;
	double delta_mom_d=0.0;
	double delta_kin_u=0.0;
	double delta_kin_d=0.0;

// kinds of parametor that I need
	int i, j, k, l,m,p;
	int t;
	int t_write=0;
	double ts = 0;
	double te = 0;

	double rx_p[NB];
	double ry_p[NB];
	double vx_p[NB];
	double vy_p[NB];
	double ax_p[NB];
	double ay_p[NB];
	double sum_vx=0.0;
	double sum_vy=0.0;

	//high speed parametor
	const double l_gx = 0.448;
	const double l_gy = l_gx;
	const int n_gx = ceil(lx / l_gx) + 6;//袖領域のため
	const int n_gy = ceil(ly / l_gy) + 6 + ceil(0.25*lx/l_gx);//袖領域のため
	const int n_g_all = n_gx * n_gy;
	int pairlist[NB][10];
	int pairlist_innit[NB][10];

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
	-3-3-3
	-2-2-2
	-1-1-1
	0 0 0   
	1 1 1 1 
	2 2 2 2 
	3 3 3 3

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

    //flywheel
    const double rf=0.25*lx;//hankei
    const int rhof=2*0.8*pow(10,pow_num);//density of fw
    const int gamma=100;//dumping by rotarion
    const double mf=M_PI*rf*rf*rhof;
    const double ia=0.5*M_PI*rhof*rf*rf*rf*rf;//inartia
	const double ria=1/ia;

    //shaft
    const double ls=pow(10,4);
    const double lsp=ls-2*rf;

	//heatwall
    double q_in=0.0;
    double q_out=0.0;
    double q_in_sum=0.0;
    double q_out_sum=0.0;

    //displacerpiston
    int mdp=2*pow(10,pow_num);
    double rmdp=1/(double) mdp;
    double dpx=0.5*lx;
    double dpy0=0.5*ly-rf;
    double dpy=dpy0;
	double dpy_b=dpy;
    double dpv =0.0;
	double dpa=0.0;

	double ddpy=0.0;
	double drdx;
	double drdx_b;
	double drdy;
	double drdy_b;

    //powerpiston
    int mpp=2*pow(10,pow_num);
    double rmpp=1/(double) mpp;
    double ppx=0.5*lx;
    double ppy0=ly;
	double ppy_max = ppy0+0.25*lx;
    double ppy=ppy0;
	double ppy_b=ppy;
    double ppv =0.0;
    double ppa=0.0;

	double dppy=0.0;
	double prdx;
	double prdx_b;
	double prdy;
	double prdy_b;

    const double mu1=ia*mdp/(ia+mdp);//mass of calculate
    const double mu2=ia*mpp/(ia+mpp);

    //connect dis-p to shaft
    double lscos=ls*(sqrt(ls*ls-rf*rf)/ls);
    //phai
    double dwx=2*rf;
	double dwx_b;
    double dwy=dpy+ls;
    double dwy_b;

    //connect pow-p to shaft
    double pwx=3*rf;
	double pwx_b;
    double pwy=lsp+ppy0;
	double pwy_b;

    //initial arg,vel,acc
    double theta=0.0;

    double omega_0=-0.01;
    double omega=omega_0;

    double alpha=0.0;

	double dtheta=0.0;
    double wy_b=lscos+dpy0;

	//work and thermalefficiency
	double w=0.0;
	double e=0.0;


	//hit or through(piston)
	int hit_piston=0;
	int through_piston=0;

	//tempreture and pressure
	double temp_d=0.0;
	double temp_u=0.0;
	double temp_all=0.0;
	double temp_lis[NWRITE];
	double press_d=0.0;
	double press_u=0.0;
	double press_all=0.0;
	double press_lis[NWRITE];
	int N_D=0;
	int N_U=0;
	int N_D_b=0;
	int N_U_b=0;

    int N_U_l_b;
    int N_U_init_sum=0;
    int N_U_l_sum=0;

	//retraint dynamics parametor
	double lambda1;
	double lambda2;
	double dpusai1[2];
	double dpusai2[2];
	double Rx1;
	double Ry1;
	double Rx2;
	double Ry2;
	double lambda1_af;
	double lambda2_af;
	double dpusai1_af[2];
	double dpusai2_af[2];
	double Rx1_af;
	double Ry1_af;
	double Rx2_af;
	double Ry2_af;

	//force for each piston
	double fp_sum=0.0;
	double fd_sum=0.0;
	double fpp=0.0;
	double fdp=0.0;

	//debug pameter
	double hss=0.0;
	double dw=0.0;
    double h1_d=0.0;

    //prepare read file
/*
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
*/

//big array

	double e_lis[NWRITE];
    double omega_lis[NWRITE];
    double theta_lis[NWRITE];
	double t_lis[NWRITE];


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


//11/30追加コンテンツ
	const double ack_rot=M_PI*0.5*0.125;
	double rot=0.0;
	int rot_num=0;
    int b_num=0;

   //12/24 DPにあつみを持たせる準備
    double dp_ty=dpy0;
	double dpy_t=dpy0+SIZE_DP*ly;
    double dpy_u=dpy0-SIZE_DP*ly;
	double dpy_u_b;
	double dpy_t_b;
    double f_left=0.05*lx;
    double f_right=0.95*lx;
    int count_fzone=0;
    int count_buffer[NB];

    for (i=0;i<NB;i++){
        count_buffer[i]=-1;
    }

    for (i=0;i<NB;i++){
        if(f_left-0.5 > rx_p[i] || f_right+0.5 < rx_p[i] ||dpy_t+0.5 < ry_p[i] || dpy_u-0.5 > ry_p[i]){
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
    double *vx_b;
    double *vy_b;
    double *ax_b;
    double *ay_b;

    int xxx,xx=0;
    rx=(double*)malloc(sizeof(double)*count_fzone);
    ry=(double*)malloc(sizeof(double)*count_fzone);
    vx=(double*)malloc(sizeof(double)*count_fzone);
    vy=(double*)malloc(sizeof(double)*count_fzone);
    ax=(double*)malloc(sizeof(double)*count_fzone);
    ay=(double*)malloc(sizeof(double)*count_fzone);

    rx_b=(double*)malloc(sizeof(double)*count_fzone);
    ry_b=(double*)malloc(sizeof(double)*count_fzone);
    vx_b=(double*)malloc(sizeof(double)*count_fzone);
    vy_b=(double*)malloc(sizeof(double)*count_fzone);
    ax_b=(double*)malloc(sizeof(double)*count_fzone);
    ay_b=(double*)malloc(sizeof(double)*count_fzone);

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
    /////////////////////////////
    int N_U_init[N];
    int N_U_l[N];
	int N_U_l_b_lis[N];
/*
	double temp_d_lis[NWRITE];
	double temp_u_lis[NWRITE];
	double press_d_lis[NWRITE];
	double press_u_lis[NWRITE];
    int ND_lis[NWRITE];
	int NU_lis[NWRITE];
    double ppy_lis[NWRITE];
    double dpy_lis[NWRITE];
	double total_kin_lis[NWRITE];
	double pot_lis[NWRITE];
	double total_e_lis[NWRITE];
*/

	double temp_part_lis_up[NWRITE][PARTNUM/2];
	double temp_part_lis_down[NWRITE][PARTNUM/2];
	double press_part_lis_up[NWRITE][PARTNUM/2];
	double press_part_lis_down[NWRITE][PARTNUM/2];
	for (i = 0;i < N;i++) {
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


///read vx,vy,file
/*
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
			ax[k] = 0.0;
			ay[k] = 0.0;
			k = k + 1;
		}
	}
*/
    //initialize N_list
    for(i=0;i<N;i++){
        N_U_init[i]=-1;
        N_U_l[i]=-1;
    }

	for (i = 0;i < N;i++) {
		kin0 = kin0 + (vx[i] * vx[i] + vy[i] * vy[i]);
	}
	kin0 = kin0 * 0.5;

	gmap_create(N, rx, ry, l_gx, l_gy, n_gx, n_gy, neighbor_list_row, neighbor_list_col, neighbor_len, pairlist, lx, g_map);
	l=0;
	m=0;
	for (i=3;i<n_gy-3;i++){
		for(j=3;j<n_gx-3;j++){
			if(g_map[i][j]!=-1){
				l+=1;
			}
		}
	}

	for(i<0;i<NWRITE;i++){
		for(j=0;j<PARTNUM/2;j++){
			temp_part_lis_up[i][j]=0.0;
			press_part_lis_down[i][j]=0.0;
			temp_part_lis_down[i][j]=0.0;
			press_part_lis_up[i][j]=0.0;
		}
	}

	printf("particlenum:%d\n",l);
	printf("kin0:%lf\n",kin0);

ppy_b=ppy0;

velocity(N,ry,vx,vy,dpy,&temp_d,&temp_u,&N_D_b,&N_U_b,N_U_init,&temp_all);
velocity(N,ry,vx,vy,dpy,&temp_d,&temp_u,&N_D,&N_U,N_U_l,&temp_all);
N_U_init_sum=int_array_sum(N,N_U_init);


//animationfile open
rx_list=fopen("./plotdata/anidat/makefig/rx_list.dat","w");
ry_list=fopen("./plotdata/anidat/makefig/ry_list.dat","w");
dwx_file=fopen("./plotdata/anidat/makefig/dwx.dat","w");
dwy_file=fopen("./plotdata/anidat/makefig/dwy.dat","w");
pwx_file=fopen("./plotdata/anidat/makefig/pwx.dat","w");
pwy_file=fopen("./plotdata/anidat/makefig/pwy.dat","w");
dp=fopen("./plotdata/anidat/makefig/dp.dat","w");
py=fopen("./plotdata/anidat/makefig/pp.dat","w");
theta_file=fopen("./plotdata/anidat/makefig/theta.dat","w");


rf_file=fopen("./plotdata/anidat/makefig/rf.dat","w");
fprintf(rf_file,"%lf",rf);
fclose(rf_file);
dpx_file=fopen("./plotdata/anidat/makefig/dpx.dat","w");
fprintf(dpx_file,"%lf",dpx);
fclose(dpx_file);

sprintf(name_kin,"%s%d%s",plotdata,(int) temp_l,c_kin);
sprintf(name_pot,"%s%d%s",plotdata,(int) temp_l,c_pot);
sprintf(name_tot,"%s%d%s",plotdata,(int) temp_l,c_tot);
sprintf(name_ene,"%s%d%s",plotdata,(int) temp_l,c_ene);
sprintf(name_dis,"%s%d%s",plotdata,(int) temp_l,c_dis);
sprintf(name_pow,"%s%d%s",plotdata,(int) temp_l,c_pow);
sprintf(name_ome,"%s%d%s",plotdata,(int) temp_l,c_ome);
sprintf(name_the,"%s%d%s",plotdata,(int) temp_l,c_the);
sprintf(name_e,"%s%d%s",plotdata,(int) temp_l,c_e);
sprintf(name_press_d,"%s%d%s",plotdata,(int) temp_l,c_press_d);
sprintf(name_press_u,"%s%d%s",plotdata,(int) temp_l,c_press_u);
sprintf(name_press,"%s%d%s",plotdata,(int) temp_l,c_press);
sprintf(name_temp,"%s%d%s",plotdata,(int) temp_l,c_temp);
sprintf(name_press_part,"%s%d%s",plotdata,(int) temp_l,c_press_part);
sprintf(name_temp_part,"%s%d%s",plotdata,(int) temp_l,c_temp_part);
sprintf(name_all_condition,"%s%d%s",plotdata,(int) temp_l,c_all_condition);
all_condition_file=fopen(name_all_condition,"w");
//-------------start mainroop------------------
	ts = omp_get_wtime();
	for (t = 1;t <= NSTEPS;t++)
	{
		//save before data
		for(i=0;i<N;i++){
			rx_b[i]=rx[i];
			ry_b[i]=ry[i];
			vy_b[i]=vy[i];
            ax_b[i]=ax[i];
			ay_b[i]=ay[i];
			N_U_l_b_lis[i]=N_U_l[i];
		}





			dpy_b=dpy;	
			dpy_u_b=dpy_u;
			dpy_t_b=dpy_t;
			dwx_b=dwx;
			dwy_b=dwy;
			ppy_b=ppy;	
			N_U_l_b=N_U_l_sum;

		//debug

		//initialize
		total_kin = 0.0;
		total_e = 0.0;
		if ((t-1)%NDEVIDE == 0){
			t_write+=1;
			t_lis[t_write-1] = t * h;
			//圧力などの定義
			partition(ry,vx,vy,ppy,dpy,temp_part_lis_down,temp_part_lis_up,press_part_lis_down,press_part_lis_up,t_write,N);
		}

		//for animaiton

            if (t%NDEVIDE == 0){
                for (i=0;i<N;i++){
                    fprintf(rx_list,"%lf    ",rx[i]);
                }
                fprintf(rx_list,"\n");

                for (i=0;i<N;i++){
                    fprintf(ry_list,"%lf    ",ry[i]);
                }
                fprintf(ry_list,"\n");
				fprintf(dp,"%lf\n",dpy);
				fprintf(py,"%lf\n",ppy);
                fprintf(pwx_file,"%lf\n",pwx);
                fprintf(pwy_file,"%lf\n",pwy);
                fprintf(dwx_file,"%lf\n",dwx);
                fprintf(dwy_file,"%lf\n",dwy);
                fprintf(theta_file,"%lf\n",theta);               
            }

		//verlet
		for (i = 0;i < N;i++) {
			rx[i] = rx[i] + vx[i] * h + ax[i] * h2;
			ry[i] = ry[i] + vy[i] * h + ay[i] * h2;
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}

		//make gmap and pairlist
		for (i = 0;i < n_gy;i++) {
			for (j = 0;j < n_gx;j++) {
				g_map[i][j] = -1;
			}
		}

		for (i = 0;i < N;i++) {
			for (j = 0;j < 10;j++) {
				pairlist[i][j] = -1;
			}
		}



	if(judge_max(ry,N,&b_num)>ppy || judge_max(rx,N,&b_num)>2*lx){
		printf("problem over p num = %d\nthis value has rx:%lf,ry:%lf,ax:%lf,ay:%lf,\nand t is %d",b_num,rx[b_num],ry[b_num],ax_b[b_num],ay_b[b_num],t);
		break;
	}

	if(judge_mini(ry,N,&b_num)<0.0 || judge_mini(rx,N,&b_num)<-2*lx){
		printf("problem under p num = %d\nthis value has rx:%lf,ry:%lf,ax:%lf,ay:%lf,\nand t is %d",b_num,rx[b_num],ry[b_num],ax_b[b_num],ay[b_num],t);
		break;
	}

		gmap_create(N, rx, ry, l_gx, l_gy, n_gx, n_gy, neighbor_list_row, neighbor_list_col, neighbor_len, pairlist, lx, g_map);

		//cliculate force		
		force(N, rx, ry, ax, ay, lx, ly, pairlist, rc2, &pot, &pot_ij, &mini);

		//second verlet
		for (i = 0;i < N;i++) {
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}

		sum_vx=double_array_sum(N,vx);
		sum_vy=double_array_sum(N,vy);

		//verocity verlet for piston
		ddpy=sqrt(pow((dwx-dpx),2)+pow((dwy-dpy),2))-ls;
		dppy=sqrt(pow((pwx-ppx),2)+pow((pwy-ppy),2))-lsp;
		dwx_b=dwx;
	    dwy_b=dwy;
		pwx_b=pwx;
		pwy_b=pwy;
		dpy+=dpv*h+0.5*dpa*h*h;
		dpy_u=dpy-SIZE_DP*ly;
		dpy_t=dpy+SIZE_DP*ly;
		ppy+=ppv*h+0.5*ppa*h*h;
		theta+=omega*h+0.5*alpha*h*h;
		dtheta=omega*h+0.5*alpha*h*h;

		//restraint dynamics

		//about disp-p
		dwx=2*rf+rf*sin(theta);
		dwy=ls+dpy0+rf*(1-cos(theta));
		drdx=dwx-dpx;
		drdx_b=dwx_b-dpx;
		drdy=dwy-dpy;
		drdy_b=dwy_b-dpy_b;
		
		//about power-p
		pwx=2*rf+rf*cos(theta);
		pwy=lsp+ppy0+rf*sin(theta);
		prdx=pwx-ppx;
		prdx_b=pwx_b-ppx;
		prdy=pwy-ppy;
		prdy_b=pwy_b-ppy_b;

		//restraint dynamics(main)
		lambda1=-1e6*mu1*(ls*ls-(drdx*drdx+drdy*drdy))/(4*(drdx*drdx_b+drdy*drdy_b));
		lambda2=-1e6*mu2*(lsp*lsp-(prdx*prdx+prdy*prdy))/(4*(prdx*prdx_b+prdy*prdy_b));
		dpusai1[0]=-2*(dwx_b-dpx);
		dpusai1[1]=-2*(dwy_b-dpy_b);
		dpusai2[0]=-2*(pwx_b-ppx);
		dpusai2[1]=-2*(pwy_b-ppy_b);
		Rx1=-1*lambda1*dpusai1[0];//Rx is forced from disp-p
		Ry1=-1*lambda1*dpusai1[1];
		Rx2=-1*lambda2*dpusai2[0];//Ry is forced from power-p
		Ry2=-1*lambda2*dpusai2[1];
		dpy+=0.5*h*h*Ry1*rmdp;//include restraint force
		dpy_u=dpy-SIZE_DP*ly;
		dpy_t=dpy+SIZE_DP*ly;
		ppy+=0.5*h*h*Ry2*rmpp;
		theta+=0.5*h*h*(-Ry2*cos(theta))*rf*ria;

		//update
		dwx=2*rf+rf*sin(theta);
		dwy=ls+dpy0+rf*(1-cos(theta));
		pwx=2*rf+rf*cos(theta);
		pwy=lsp+ppy0+rf*sin(theta);
		dpv+=dpa*h*0.5+0.5*h*Ry1*rmdp;
		ppv+=ppa*h*0.5+0.5*h*Ry2*rmpp;
		omega+=0.5*h*alpha+0.5*(-Ry2*cos(theta)*h*rf)*ria;

		//piston force
		heatwall(h,N,ry,ry_b,vy,&q_in,&q_out,ppy,ppv,temp_l,temp_h,&fpp,ly,&hss,&dw);
		piston_move_u(N,ry,ry_b,vy,vy_b,ay_b,h,h_rev,&q_in,&q_out,q_in_sum,q_out_sum,dpy_t,dpy_t_b,&dpv,&hit_piston,&through_piston,&fdp,temp_l,temp_h,&h1_d,momentum_u,kin_u,&up_hit,&up_through,mdp,&delta_mom_d,&delta_kin_d,f_right,f_left,rx);
		piston_move_d(N,ry,ry_b,vy,vy_b,ay_b,h,h_rev,&q_in,&q_out,q_in_sum,q_out_sum,dpy_u,dpy_u_b,&dpv,&hit_piston,&through_piston,&fdp,temp_l,temp_h,&h1_d,momentum_d,kin_d,&down_hit,&down_through,mdp,&delta_mom_u,&delta_kin_u,f_right,f_left,rx);
//		fix(rx,ry,vx,vy,f_left,f_right,dpy_t,dpy_u,N);
		//It can be better
		ppa=fpp*rmpp;
		dpa=fdp*rmdp;
		dpv+=0.5*h*dpa;
		ppv+=0.5*h*ppa;
		alpha=-gamma*omega*ria;
		omega+=0.5*h*alpha;
		lambda1_af=1e3*(dpv-rf*omega*sin(theta))*(dpy-dwy)*mu1/(ls*ls);
    	lambda2_af=1e3*(ppv-rf*omega*cos(theta))*(ppy-pwy)*mu2/(lsp*lsp);
		dpusai1_af[0]=-2*(dwx-dpx);
		dpusai1_af[1]=-2*(dwy-dpy);
		dpusai2_af[0]=-2*(pwx-ppx);
		dpusai2_af[1]=-2*(pwy-ppy);
		Rx1_af=-1*lambda1_af*dpusai1_af[0];
		Ry1_af=-1*lambda1_af*dpusai1_af[1];
		Rx2_af=-1*lambda2_af*dpusai2_af[0];
		Ry2_af=-1*lambda2_af*dpusai2_af[1];
		dpv+=0.5*h*(Ry1_af)*rmdp;
		ppv+=0.5*h*(Ry2_af)*rmpp;
		omega+=0.5*h*(-Ry2_af*cos(theta))*rf*ria;

		//boundary condition 
		boundary(N, rx, lx);
		forbid_area(rx,ry,rx_b,vx,vy,f_left,f_right,dpy_t,dpy_u,N);
///////temprature and pressure onthepiston under the piston/////////

		velocity(N,ry,vx,vy,dpy,&temp_d,&temp_u,&N_D,&N_U,N_U_l,&temp_all);

		press_all=temp_all*N/(lx*ppy);
		press_d=temp_d*N_D/(lx*dpy);
		press_u=temp_u*N_U/(lx*(ppy-dpy));
        N_U_l_sum=int_array_sum(N,N_U_l);


		//caliclate kin-energy
		for (i = 0;i < N;i++) {
			total_kin += (vx[i] * vx[i] + vy[i] * vy[i]);
		}
		total_kin = total_kin * 0.5;
		total_e = total_kin + pot;

			if (t % nout == 0){
				printf("Time:%lf,Kinetic Energy:%lf,Potential Energy:%lf,Total Energy:%lf,momentum_vx:%lf,momentum_vy:%lf\n", t * h, total_kin, pot, total_e,sum_vx,sum_vy);
	}
		//work and thermal efficiency
		w+=gamma*rf*omega*dtheta;
		e=w/(0.5*ia*omega_0*omega_0+q_in);

		rot=fabs(theta)-rot_num*ack_rot;
		if (rot>ack_rot){
			rot-=fabs(ack_rot);
			fprintf(all_condition_file,"%lf    %lf    %lf    %lf    %lf    %lf\n",t*h,press_all,temp_all,ppy*lx,ppy,dpy);
			rot_num+=1;
		}


		//into the big array
		if ((t-1)%NDEVIDE == 0)
		{
			e_lis[t_write-1]=e;
			omega_lis[t_write-1]=omega;
			theta_lis[t_write-1]=theta;
			press_lis[t_write-1]=press_all;
			temp_lis[t_write-1]=temp_all;

/*
		temp_d_lis[t_wirte-1]=temp_d;
		temp_u_lis[t_wirte-1]=temp_u;
		press_d_lis[t_wirte-1]=press_d;
		press_u_lis[t_wirte-1]=press_u;
		dpy_lis[t_wirte-1]=dpy;
		ppy_lis[t_wirte-1]=ppy;
		NU_lis[t_wirte-1]=N_U;
		ND_lis[t_wirte-1]=N_D;
		total_kin_lis[t_wirte-1] = total_kin;
		pot_lis[t_wirte-1] = pot;
		total_e_lis[t_wirte-1] = total_e;
*/		
		}
/*
	if(judge_max(ry,N)>ppy){
		printf("over pp\n");
		break;
	}
*/
}
//end mainloop
te = omp_get_wtime();
printf("time cost is %lf seconds\n", te - ts);
printf("process time is:%lf\n",te-ts);
printf("hit is :%d,through is %d\n",hit_piston,through_piston);
printf("up_hit:%d,down_hit:%d\n",up_hit,down_hit);
printf("up_through:%d,down_through:%d\n",up_through,down_through);

//free gmap
for(i=0;i<n_gy;i++){
	free(g_map[i]);
}
free(g_map);


//file output


efile = fopen(name_e,"w");
omega_file=fopen(name_ome,"w");
//theta_file=fopen(name_the,"w");
press_file=fopen(name_press,"w");
temp_file=fopen(name_temp,"w");
part_press_file=fopen(name_press_part,"w");
part_temp_file=fopen(name_temp_part,"w");
/*
kin_file = fopen(name_kin,"w");
pot_file = fopen(name_pot,"w");
tot_file = fopen(name_tot,"w");
enegy_file=fopen(name_ene,"w");
dis_p = fopen(name_dis,"w");
pow_p = fopen(name_pow,"w");
press_d_file=fopen(name_press_d,"w");
press_u_file=fopen(name_press_u,"w");
*/


//write in file
for(i=0;i<NWRITE;i++){
		t_lis[i]-=h;
		fprintf(omega_file,"%lf   %lf\n",t_lis[i],omega_lis[i]);
//		fprintf(theta_file,"%lf    %lf\n",t_lis[i],theta_lis[i]);
		fprintf(efile,"%lf    %lf\n",t_lis[i],e_lis[i]);
		fprintf(part_temp_file,"%lf    ",t_lis[i]);
		fprintf(part_press_file,"%lf    ",t_lis[i]);
		for (j=0;j<PARTNUM/2;j++){
			fprintf(part_temp_file,"%lf    ",temp_part_lis_down[i][j]);
			fprintf(part_press_file,"%lf    ",press_part_lis_down[i][j]);
			fprintf(part_temp_file,"%lf    ",temp_part_lis_up[i][j]);
			fprintf(part_press_file,"%lf    ",press_part_lis_up[i][j]);
		}
		fprintf(part_temp_file,"\n");
		fprintf(part_press_file,"\n");
		fprintf(press_file,"%lf    %lf\n",t_lis[i],press_lis[i]);
		fprintf(temp_file,"%lf    %lf\n",t_lis[i],temp_lis[i]);
/*		
		fprintf(kin_file,"%lf	%lf\n",t_lis[i],total_kin_lis[i]);
		fprintf(pot_file,"%lf	%lf\n",t_lis[i],pot_lis[i]);
		fprintf(tot_file,"%lf	%lf\n",t_lis[i],total_e_lis[i]);
		fprintf(enegy_file,"%lf    %lf    %lf    %lf\n",t_lis[i],total_kin_lis[i],pot_lis[i],total_e_lis[i]);
		fprintf(dis_p,"%lf	%lf	  \n",t_lis[i],dpy_lis[i]);
		fprintf(pow_p,"%lf	%lf	  \n",t_lis[i],ppy_lis[i]);
		fprintf(press_d_file,"%lf    %lf\n",t_lis[i],press_d_lis[i]);
		fprintf(press_u_file,"%lf    %lf\n",t_lis[i],press_u_lis[i]);
*/

	
}

//file close

fclose(all_condition_file);

fclose(omega_file);
//fclose(theta_file);
fclose(efile);
fclose(press_file);
fclose(temp_file);
fclose(part_press_file);
fclose(part_temp_file);


lxly=fopen("./plotdata/anidat/makefig/lxly.dat","w");
fprintf(lxly,"%lf    %lf    %lf     %d",lx,ppy_max,temp_l,PARTNUM);
fclose(lxly);


/*
fclose(kin_file);
fclose(pot_file);
fclose(tot_file);
fclose(enegy_file);
fclose(dis_p);
fclose(pow_p);
fclose(press_d_file);
fclose(press_u_file);
*/



//close particle animation
fclose(rx_list);
fclose(ry_list);
fclose(dp);
fclose(py);
	return 0;

}






//prototype function after write
void gmap_create(int NP, double RX[], double RY[], double L_GX, double L_GY, int N_GX, int N_GY, int NEIGHBOR_ROW[], int NEIGHBOR_COL[], int NEITGHBOR_LEN, int PAIRLIST[][10], double LX, int **G_MAP) {
	int i, j, k;
	int gx_map;
	int gy_map;
	int partcle_counter;

	for (i = 0;i < NP;i++) {

		if (RX[i] >= LX) {
			RX[i] = RX[i] - LX;
		}

		else if (RX[i] <= 0.0) {
			RX[i] = RX[i] + LX;
		}

		gx_map = (int)(RX[i] / L_GX) + 3;
		gy_map = N_GY - 4 - (int)(RY[i] / L_GY);//-1+3-6 (int)(RY[i] / L_GY) + 3
		G_MAP[gy_map][gx_map] = i;

	}



	for (i = 0;i < N_GY;i++) {
		for (j = 0;j < 3;j++) {
			G_MAP[i][j] = G_MAP[i][N_GX - 6 + j];
			G_MAP[i][j + N_GX - 3] = G_MAP[i][3 + j];
		}
	}

	for (i = 3;i < N_GY -3;i++)
	 {
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
	int roop_num;
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
			}else {

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




void heatwall(double H,int NP,double RY[],double RY_B[],double VY[],double *Q_IN,double *Q_OUT ,double PPY,double PPV,double TEMP_L,double TEMP_H,double *FPP,double LY,double *h_ss,double *d_w){
	double dq=0.0;
	double fpu=0.0;
	double fpu_sum=0.0;
	double h1=0.0;
	double vy_l=0.0;
	int i;
	for(i=0;i<NP;i++){


		if(RY[i]<=0.5){
			
			*d_w=RY_B[i]-0.5;
			h1=(*d_w)/fabs(VY[i]);
			*h_ss=H-h1;
			RY[i]=RY_B[i]+VY[i]*h1;

			vy_l=VY[i];
			VY[i]=canon_b(TEMP_L);//canon_b(TEMP_L)
			RY[i]+=VY[i]*(*h_ss);
			dq=0.5*(VY[i]*VY[i]-vy_l*vy_l);

			if(dq>=0.0){
				*Q_IN+=dq;
			}else{
				*Q_OUT+=dq;
			}

		}



		if(RY[i]>=PPY-0.5){
			*d_w=(PPY-0.5)-RY_B[i];
			h1=*d_w/fabs(VY[i]-PPV);
			*h_ss=H-h1;
			RY[i]=RY_B[i]+VY[i]*h1;
			vy_l=VY[i];
			VY[i]=-canon_b(TEMP_H);
			RY[i]+=VY[i]*(*h_ss);
			fpu=-(VY[i]-vy_l)/H;
			fpu_sum+=fpu;
			dq=0.5*(VY[i]*VY[i]-vy_l*vy_l);

			if(dq>=0.0){
				*Q_IN+=dq;
			}else{
				*Q_OUT+=dq;
			}
		}

	}

		*FPP=fpu_sum;
}




void piston_move_u(int NP,double RY[],double RY_B[],double VY[],double VY_B[],double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,double *H1_D,double Momentum_u[2],double Kin_u[2],int *k,int *j,double MDP,double *delta_mom,double *delta_kin,double F_RIGHT,double F_LEFT,double RX[]){
	double dq=0.0;
	double fdu=0.0;
	double fdu_sum=0.0;
	double fdd=0.0;
	double fdd_sum=0.0;
	double h_ss=0.0;
	double d_w=0.0;
	double dd=0.0;
	double h1=0.0;
	double vy_l=0.0;
	double rand_num;
	int i;
	double re_vel=0.0;

	for(i=0;i<NP;i++){


		if((F_RIGHT+0.5)>RX[i] && (F_LEFT-0.5)<RX[i]){
			if(RY_B[i]>=DPY_B+0.5 && RY[i]<DPY+0.5){
				//devide step h -> h1,h2 h1
				rand_num=genrand_real3();
				d_w=RY_B[i]-(DPY_B+0.5);
				re_vel=(VY_B[i]-(*DPV))*(VY_B[i]-(*DPV));
				re_vel=sqrt(re_vel);
				h1=d_w/re_vel;
				if (h1<0){printf("\nup:%lf\n",h1);}
				*H1_D=h1;
				h_ss=fabs(H-h1);
				RY[i]=RY_B[i]+VY_B[i]*h1;
				vy_l=VY[i];
				//devide step h -> h1,h2 h2
			if (rand_num<PROBABIRITY){
					*HIT_PISTON+=1;
					*k+=1;
					Momentum_u[0]=VY[i]+MDP*(*DPV);
					Kin_u[0]=0.5*VY[i]*VY[i]+0.5*MDP*(*DPV)*(*DPV);

					VY[i]=-((MDP-1)/(MDP+1))*VY_B[i]+((2*MDP)/(MDP+1))*(*DPV);
				*DPV=((MDP-1)/(MDP+1))*(*DPV)+(2/(MDP+1))*vy_l;

					Momentum_u[1]=VY[i]+MDP*(*DPV);
					Kin_u[1]=0.5*VY[i]*VY[i]+0.5*MDP*(*DPV)*(*DPV);

					*delta_mom+=Momentum_u[1]-Momentum_u[0];
					*delta_kin+=Kin_u[1]-Kin_u[0];
					*delta_mom=(*delta_mom)/(*k);
					*delta_kin=(*delta_kin)/(*k);

					dd=VY[i]*h_ss;
					RY[i]=RY[i]+dd;
					fdu=-(VY[i]-VY_B[i])/H;
					fdu_sum+=fdu;
					dq=0.5*(VY[i]*VY[i]-VY_B[i]*VY_B[i]);
						if(dq>=0.0){
							*Q_IN+=dq;
						}else{
							*Q_OUT+=dq;
						}


				}else{
					*j+=1;
					*THROUGH_PISTON+=1;

					dd=VY[i]*h_ss;
					RY[i]=RY[i]+dd;
					
				}

			}

		}


		//dendfor
		}
	*FDP=fdd_sum+fdu_sum;
}
	
void piston_move_d(int NP,double RY[],double RY_B[],double VY[],double VY_B[],double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,double *H1_D,double Momentum_d[2],double Kin_d[2],int *kk,int *jj,double MDP,double *delta_mom,double *delta_kin,double F_RIGHT,double F_LEFT,double RX[]){
	double dq=0.0;
	double fdu=0.0;
	double fdu_sum=0.0;
	double fdd=0.0;
	double fdd_sum=0.0;
	double h_ss=0.0;
	double d_w=0.0;
	double dd=0.0;
	double h1=0.0;
	double vy_l=0.0;
	double z=0.0;
	int i;
	double re_vel=0.0;
    for(i=0;i<NP;i++){
		if((F_RIGHT+0.5)>RX[i] && (F_LEFT-0.5)<RX[i]){
			if(RY_B[i]<=DPY_B-0.5 && RY[i]>DPY-0.5){
				z=genrand_real3();
				d_w=(DPY_B-0.5)-RY_B[i];
				re_vel=(VY_B[i]-(*DPV))*(VY_B[i]-(*DPV));
				re_vel=sqrt(re_vel);
				h1=d_w/re_vel;
				h_ss=fabs(H-h1);

				RY[i]=RY_B[i]+VY_B[i]*h1;
				DPY=DPY_B+(*DPV)*h1;
				if (z<PROBABIRITY){
					*kk+=1;
					*HIT_PISTON+=1;
					vy_l=VY[i];
					Momentum_d[0]=VY[i]+MDP*(*DPV);
					Kin_d[0]=0.5*VY[i]*VY[i]+0.5*MDP*(*DPV)*(*DPV);

					VY[i]=-((MDP-1)/(MDP+1))*VY_B[i]+((2*MDP)/(MDP+1))*(*DPV);
					*DPV=((MDP-1)/(MDP+1))*(*DPV)+(2/(MDP+1))*vy_l;

					Momentum_d[1]=VY[i]+MDP*(*DPV);
					Kin_d[1]=0.5*VY[i]*VY[i]+0.5*MDP*(*DPV)*(*DPV);

					*delta_mom+=Momentum_d[1]-Momentum_d[0];
					*delta_kin+=Kin_d[1]-Kin_d[0];
					*delta_mom=(*delta_mom)/(*kk);
					*delta_kin=(*delta_kin)/(*kk);

					dd=VY[i]*h_ss;
					if (VY[i]>0){
					}
					else{RY[i]=RY[i]+dd;}
					if((*DPV)<0.0 && RY[i]>DPY-0.5){
						printf("h_ss :%lf,d_w :%lf,re_vel:%lf,h1:%lf,vy_b:%lf,dpv:%lf\n",h_ss,d_w,re_vel,h1,VY_B[i],*DPV);
						}
					fdd=-(VY[i]-vy_l)*H_REV;
					fdd_sum+=fdu;
					dq=0.5*(VY[i]*VY[i]-vy_l*vy_l);
					if(dq>=0.0){
						*Q_IN+=dq;
					}else{
						*Q_OUT+=dq;
					}


				}else{
					*jj+=1;
					*THROUGH_PISTON+=1;
					dd=VY[i]*h_ss;
					RY[i]=RY[i]+dd;                  
				}

			}

		}




    }
	*FDP=fdd_sum+fdu_sum;

}

void velocity(int NP,double RY[],double VX[],double VY[],double DPY,double *TEMP_D,double *TEMP_U,int *ND,int*NU ,int N_U_list[],double *TEMP){
	*TEMP_D=0.0;
	*TEMP_U=0.0;
	*TEMP=0.0;
	*NU=0.0;
	*ND=0.0;
	int i;
    int k=0;
    int *p;
    p=N_U_list;
	for(i=0;i<NP;i++){
		*TEMP+=0.5*(VX[i]*VX[i]+VY[i]*VY[i]);
		if(RY[i]<DPY-0.5){
			*TEMP_D+=0.5*(VX[i]*VX[i]+VY[i]*VY[i]);
			*ND+=1;
		}else{
			*TEMP_U+=0.5*(VX[i]*VX[i]+VY[i]*VY[i]);
			*NU+=1;
			N_U_list[i]=i;
		}
	}
	*TEMP=*TEMP/NP;
	*TEMP_D=*TEMP_D/(*ND);
	*TEMP_U=*TEMP_U/(*NU);
}


void partition(double RY[],double VX[] , double VY[] ,double PPY,double DPY,double TEMP_PART_LIST_DOWN[NWRITE][PARTNUM/2],double TEMP_PART_LIST_UP[NWRITE][PARTNUM/2],double PRESS_PART_LIST_DOWN[NWRITE][PARTNUM/2],double PRESS_PART_LIST_UP[NWRITE][PARTNUM/2],int T_WRITE,int NP){
	int i,j,k,l;
	double temp;
	double partnum= PARTNUM/2;
 	double re_partnum=(double) 1/partnum;
	double part,part_plus;
	double up_pow_ppy=((PPY-DPY)*re_partnum)*((PPY-DPY)*re_partnum);
	double down_pow_ppy=(DPY*re_partnum)*(DPY*re_partnum);
	double up_y;
	double up_part,down_part,up_part_plus,down_part_plus;
	double up_block,down_block;
	double temp_up,temp_down;
	for(i=0;i<partnum;i++){
		k=0;
		l=0;
		up_y=PPY-DPY;
		up_block=up_y*re_partnum;
		down_block=DPY*re_partnum;
		up_part=((double) (i/partnum))*up_y+DPY;
		up_part_plus=((double) (i+1)/partnum)*up_y+DPY;
		down_part=((double) (i/partnum))*DPY;
		down_part_plus=((double) (i+1)/partnum)*DPY;
		for(j=0;j<NP;j++){
			if(up_part<=RY[j] && up_part_plus>RY[j]){
				k+=1;
				TEMP_PART_LIST_UP[T_WRITE-1][i]+=0.5*(VX[j]*VX[j]+VY[j]*VY[j]);
				}
			if(down_part<=RY[j] && down_part_plus>RY[j]){
				l+=1;
				TEMP_PART_LIST_DOWN[T_WRITE-1][i]+=0.5*(VX[j]*VX[j]+VY[j]*VY[j]);
				}
			}
		TEMP_PART_LIST_UP[T_WRITE-1][i]=TEMP_PART_LIST_UP[T_WRITE-1][i]/k;
		TEMP_PART_LIST_DOWN[T_WRITE-1][i]=TEMP_PART_LIST_DOWN[T_WRITE-1][i]/l;
		temp_up=TEMP_PART_LIST_UP[T_WRITE-1][i];
		temp_down=TEMP_PART_LIST_DOWN[T_WRITE-1][i];
		PRESS_PART_LIST_DOWN[T_WRITE-1][i]=(l/down_pow_ppy)*temp_down;
		PRESS_PART_LIST_UP[T_WRITE-1][i]=(k/up_pow_ppy)*temp_up;
	}
}