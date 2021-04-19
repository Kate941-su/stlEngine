//g_mapの配列をmallocすることでもっといいコードがかけるかもしれないがc++で実装するのもあり
//とりあえず静的な二次元配列を確保しておく

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "MT.h"
#include <time.h>

#define Nx 15
#define Ny 30
#define N (Nx*Ny)
#define NSTEPS 10000


//prototype function
void gmap_create(int NP, double RX[], double RY[], double L_GX, double L_GY, int N_GX, int N_GY, int NEIGHBOR_ROW[], int NEIGHBOR_COL[], int NEITGHBOR_LEN, int PAIRLIST[][10], double LX, int **G_MAP);

void force(int NP, double RX[], double RY[], double AX[], double AY[], double LX, double LY, int PAIRLIST[][10], double RC2, double* POT, double* POT_IJ, double* MINI);

void heatwall(double H,int NP,double RY[],double RY_B[],double VY[],double *Q_IN,double *Q_OUT ,double PPY,double PPV,double TEMP_L,double TEMP_H,double *FPP,double LY,double *h_ss,double *d_w);

void piston_move(int NP,double RY[],double RY_B[],double VY[],double VY_B[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H);

void velocity(int NP,double RY[],double VX[],double VY[],double DPY,double *TEMP_D,double *TEMP_U,int *ND,int*NU );


//return max
double judge_max(double array[],int len_array){
	int i;
	double max=0.0;
	for(i=0;i<len_array;i++){
		if (array[i]>max){
			max=array[i];
		}
	}
	return max;
}

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

	//for animation
	FILE *rx_list;
	FILE *ry_list;

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
	const double h_rev=1/h;
	const double h2 = 0.5 * h * h;
	const int nout = NSTEPS / 10;
	const int t_end = NSTEPS * h;

	double kin0=0;
	double pot = 0;
	double pot_ij = 0;
	double pot_lis[NSTEPS];
	double total_kin=0;
	double total_kin_lis[NSTEPS];
	double total_e=0;
	double total_e_lis[NSTEPS];


	int i, j, k, l,m,p;
	int t;
	double ts = 0;
	double te = 0;
	double t_lis[NSTEPS];
	double rx[N];
	double rx_b[N];
	double ry[N];
	double ry_b[N];
	double vx[N];
	double vy[N];
	double vy_b[N];
	double ax[N];
	double ay[N];
	double ay_b[N];


	//high speed parametor
	const double l_gx = 0.56;
	const double l_gy = l_gx;
	const int n_gx = ceil(lx / l_gx) + 6;//袖領域のため
	const int n_gy = ceil(ly / l_gy) + 6 + ceil(0.25*lx/l_gx);//袖領域のため
	const int n_g_all = n_gx * n_gy;
	int pairlist[N][10];

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

    //flywheel
    const double rf=0.25*lx;//hankei
    const int rhof=0.8*pow(10,4);//density of fw
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
    double q_in_lis[NSTEPS];
    double q_out_lis[NSTEPS];

    //displacerpiston
    const int mdp=pow(10,3);
    const double rmdp=1e-3;
    double dpx=0.5*lx;
    double dpy0=0.5*ly-rf;
    double dpy=dpy0;
	double dpy_b=dpy;
    double dpv =0.0;
	double dpa=0.0;
    double dpy_lis[NSTEPS];
	double ddpy=0.0;
	double drdx;
	double drdx_b;
	double drdy;
	double drdy_b;

    //powerpiston
    const int mpp=pow(10,3);
    const double rmpp=1e-3;
    double ppx=0.5*lx;
    double ppy0=ly;
	double ppy_max = ppy0+0.25*lx;
    double ppy=ppy0;
	double ppy_b=ppy;
    double ppv =0.0;
    double ppa=0.0;
    double ppy_lis[NSTEPS];
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
    double theta_lis[NSTEPS];
    double omega_0=-0.05;
    double omega=omega_0;
    double omega_lis[NSTEPS];
    double alpha=0.0;
    double alpha_lis[NSTEPS];
	double dtheta=0.0;
    double wy_b=lscos+dpy0;

	//work and thermalefficiency
	double w=0.0;
	double e=0.0;
	double e_lis[NSTEPS];

	//hit or through(piston)
	int hit_piston=0;
	int through_piston=0;


	//input and output thermal


	//tempreture and pressure
	double temp_d=0.0;
	double temp_u=0.0;
	int N_D=0;
	int N_U=0;
	int N_D_b=0;
	int N_U_b=0;
	double temp_d_lis[NSTEPS];
	double temp_u_lis[NSTEPS];
    int ND_lis[NSTEPS];
	int NU_lis[NSTEPS];

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

	k = 0;
	for (i = 0;i < Nx;i++) {
		for (j = 0;j < Ny;j++) {
			rx[k] = i * dx + dx * 0.5;
			ry[k] = j * dy + dy * 0.5;
			vx[k] = sqrt(temp) * normal();
			vy[k] = sqrt(temp) * normal();
			ax[k] = 0.0;
			ay[k] = 0.0;
			k = k + 1;
		}
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
	printf("particlenum:%d\n",l);

//	printf("kin0:%lf\n", kin0);
//	printf("counter:%d\n", l);
	printf("kin0:%lf\n",kin0);

ppy_b=ppy0;
dpy_b=dpy0;

velocity(N,ry,vx,vy,dpy,&temp_d,&temp_u,&N_D_b,&N_U_b);
            
rx_list=fopen("./plotdata/rx_list.dat","w");
ry_list=fopen("./plotdata/ry_list.dat","w");
dp=fopen("./plotdata/dp.dat","w");
py=fopen("./plotdata/pp.dat","w");

//-------------start mainroop------------------
	ts = omp_get_wtime();
	for (t = 1;t <= NSTEPS;t++)
	{
		//save before data
		for(i=0;i<N;i++){
			rx_b[i]=rx[i];
			ry_b[i]=ry[i];
			vy_b[i]=vy[i];
			ay_b[i]=ay[i];	
		}

			dpy_b=dpy;	
			dwx_b=dwx;
			dwy_b=dwy;
			ppy_b=ppy;	


		//debug

		//initialize
		total_kin = 0.0;
		total_e = 0.0;
		t_lis[t-1] = t * h;


		//for animaiton
        if (t>=7000){
            if (t%20 == 0){
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

            }
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


		gmap_create(N, rx, ry, l_gx, l_gy, n_gx, n_gy, neighbor_list_row, neighbor_list_col, neighbor_len, pairlist, lx, g_map);

////////debug zone///////////
		l=0;
		m=0;
		for (i=0;i<n_gy;i++)//extend n_gy-3 -> n_gy
		{
			for(j=3;j<n_gx-3;j++){
				if(g_map[i][j]!=-1){
					m+=1;
					l+=g_map[i][j];
				}
			}
		}
		p=0;
		for(i=0;i<N;i++){
			p+=i;
		}

		if(m!= N){
			printf("particlenum:%d\n",m);
			printf("kasanari\n");
			printf("fade num is %d\n",p-l);
// 			break;
		}





		//cliculate force		
		force(N, rx, ry, ax, ay, lx, ly, pairlist, rc2, &pot, &pot_ij, &mini);

		//second verlet
		for (i = 0;i < N;i++) {
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}

		//verocity verlet for piston
		ddpy=sqrt(pow((dwx-dpx),2)+pow((dwy-dpy),2))-ls;
		dppy=sqrt(pow((pwx-ppx),2)+pow((pwy-ppy),2))-lsp;
		dwx_b=dwx;
	    dwy_b=dwy;
		pwx_b=pwx;
		pwy_b=pwy;
		dpy+=dpv*h+0.5*dpa*h*h;
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
		piston_move(N,ry,ry_b,vy,vy_b,h,h_rev,&q_in,&q_out,q_in_sum,q_out_sum,dpy,dpy_b,dpv,&hit_piston,&through_piston,&fdp,temp_l,temp_h);


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


///////temprature and pressure onthepiston under the piston/////////

		velocity(N,ry,vx,vy,dpy,&temp_d,&temp_u,&N_D,&N_U);


		//caliclate kin-energy
		for (i = 0;i < N;i++) {
			total_kin += (vx[i] * vx[i] + vy[i] * vy[i]);
		}
		total_kin = total_kin * 0.5;
		total_kin_lis[t-1] = total_kin;
		pot_lis[t-1] = pot;
		total_e = total_kin + pot;
		total_e_lis[t-1] = total_e;

			if (t % nout == 0){
				printf("Time:%lf,Kinetic Energy:%lf,Potential Energy:%lf,Total Energy:%lf\n", t * h, total_kin, pot, total_e);
				printf("particlenum:%d\n",m);

			}
		//work and thermal efficiency
		w+=gamma*rf*omega*dtheta;
		e=w/(0.5*ia*omega_0*omega_0+q_in_sum);
		e_lis[t-1]=e;

		//into the array
		temp_d_lis[t-1]=temp_d;
		temp_u_lis[t-1]=temp_u;
		dpy_lis[t-1]=dpy;
		ppy_lis[t-1]=ppy;
		omega_lis[t-1]=omega;
		theta_lis[t-1]=theta;
		NU_lis[t-1]=N_U;
		ND_lis[t-1]=N_D;

	if(judge_max(ry,N)>ppy){
		printf("over pp\n");
		break;
	}

	if(N_U_b!=N_U){
		printf("Through\n");
		break;
	}	

}
te = omp_get_wtime();
printf("max ay is :%lf\n",judge_max(ay,N));
printf("max ry is :%lf\n",judge_max(ry,N));
printf("cycle num is:%d\n",t-1);

//free gmap
for(i=0;i<n_gy;i++){
	free(g_map[i]);
}
free(g_map);


//end mainloop

//file output
/*
kin_file = fopen("./plotdata/kin.dat","w");
pot_file = fopen("./plotdata/pot.dat","w");
tot_file = fopen("./plotdata/tot.dat","w");
enegy_file=fopen("./plotdata/energy.dat","w");
dis_p = fopen("./plotdata/dis_p.dat","w");
pow_p = fopen("./plotdata/pow_p.dat","w");
omega_file=fopen("./plotdata/omega.dat","w");
for(i=0;i<NSTEPS;i++){
	fprintf(kin_file,"%lf	%lf\n",t_lis[i],total_kin_lis[i]);
	fprintf(pot_file,"%lf	%lf\n",t_lis[i],pot_lis[i]);
	fprintf(tot_file,"%lf	%lf\n",t_lis[i],total_e_lis[i]);
	fprintf(enegy_file,"%lf    %lf    %lf    %lf\n",t_lis[i],total_kin_lis[i],pot_lis[i],total_e_lis[i]);
	fprintf(dis_p,"%lf	%lf	  \n",t_lis[i],dpy_lis[i]);
	fprintf(pow_p,"%lf	%lf	  \n",t_lis[i],ppy_lis[i]);
	fprintf(omega_file,"%lf   %lf\n",t_lis[i],omega_lis[i]);

}
fclose(kin_file);
fclose(pot_file);
fclose(tot_file);
fclose(enegy_file);
fclose(dis_p);
fclose(pow_p);
fclose(omega_file);

printf("time cost is %lf seconds\n", te - ts);

*/

//debugfile
rx_b_file=fopen("./plotdata/rx_b.dat","w");
ry_b_file=fopen("./plotdata/ry_b.dat","w");
rx_file=fopen("./plotdata/rx.dat","w");
ry_file=fopen("./plotdata/ry.dat","w");
lxly=fopen("./plotdata/lxly.dat","w");
dp_b=fopen("./plotdata/dp_b.dat","w");
for (i=0;i<N;i++){
	fprintf(rx_b_file,"%lf    ",rx_b[i]);
	fprintf(ry_b_file,"%lf    ",ry_b[i]);
	fprintf(rx_file,"%lf    ",rx[i]);
	fprintf(ry_file,"%lf    ",ry[i]);
}

fprintf(lxly,"%lf    %lf",lx,ppy_max);

fprintf(dp_b,"%lf",dpy_b);
fclose(rx_b_file);
fclose(ry_b_file);
fclose(rx_file);
fclose(ry_file);
fclose(lxly);
fclose(dp);
fclose(dp_b);
fclose(py);
//close particle animation
fclose(rx_list);
fclose(ry_list);


//check the gmap
/*
l = 0;
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


	printf("counter:%d\n", l);
	printf("n_gx:%d\n", n_gx);
	printf("n_gy:%d\n", n_gy);
*/

	return 0;

}






//prototype function after write
void gmap_create(int NP, double RX[], double RY[], double L_GX, double L_GY, int N_GX, int N_GY, int NEIGHBOR_ROW[], int NEIGHBOR_COL[], int NEITGHBOR_LEN, int PAIRLIST[][10], double LX, int **G_MAP) {
	int i, j, k;
	int gx_map;
	int gy_map;
	int partcle_counter;

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



//		if (RY[i] > LY - 0.5) {
//			VY[i] = -VY[i];
//		}
//		if (RY[i]<0.5){
//			VY[i]=-VY[i];
//		}


	}

		*FPP=fpu_sum;

	
}




void piston_move(int NP,double RY[],double RY_B[],double VY[],double VY_B[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H){
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


	for(i=0;i<NP;i++){
		if(RY_B[i]>DPY_B+0.5 && RY[i]<DPY+0.5){
			z=genrand_real3();
			d_w=RY_B[i]-(DPY+0.5);
			h1=d_w/fabs(VY[i]-DPV);
			h_ss=H-h1;
			RY[i]=RY_B[i]+VY_B[i]*h1;
			DPY=DPY_B+DPV*h1;
			if(0<=z<1.1){
				*HIT_PISTON+=1;
				VY[i]=-VY_B[i]+2*DPV;
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
				*THROUGH_PISTON+=1;
				VY[i]=-canon_b(TEMP_L);
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
			}
		}
		else if(RY_B[i]<DPY_B-0.5 && RY[i]>DPY-0.5){
			z=genrand_real3();
			d_w=(DPY-0.5)-RY_B[i];
			h1=d_w/fabs(VY[i]-DPV);
			h_ss=H-h1;
			RY[i]=RY_B[i]+VY_B[i]*h1;
			DPY=DPY_B+DPV*h1;
			if(0<=z<1.1){
				*HIT_PISTON+=1;
				VY[i]=-VY_B[i]+2*DPV;
				dd=VY[i]*h_ss;
				RY[i]=RY[i]+dd;
				fdd=-(VY[i]-VY_B[i])*H_REV;
				fdd_sum+=fdu;
				dq=0.5*(VY[i]*VY[i]-VY_B[i]*VY_B[i]);
				if(dq>=0.0){
					*Q_IN+=dq;
				}else{
					*Q_OUT+=dq;
				}
			}else{
				*THROUGH_PISTON+=1;
				VY[i]=canon_b(TEMP_H);
				dd=VY[i]*h_ss;
				RY[i]=RY[i]+dd;
				fdu=-(VY[i]-VY_B[i])*H_REV;
				fdu_sum+=fdu;
				dq=0.5*(VY[i]*VY[i]-VY_B[i]*VY_B[i]);
				if(dq>=0.0){
					*Q_IN+=dq;
				}else{
					*Q_OUT+=dq;
				}
			}
		}
		

	}
	*FDP=fdd_sum+fdu_sum;
}
	








void velocity(int NP,double RY[],double VX[],double VY[],double DPY,double *TEMP_D,double *TEMP_U,int *ND,int*NU ){
	*TEMP_D=0.0;
	*TEMP_U=0.0;
	*NU=0.0;
	*ND=0.0;
	int i;
	for(i=0;i<NP;i++){
		if(RY[i]<DPY){
			*TEMP_D+=0.5*(VX[i]*VX[i]+VY[i]*VY[i]);
			*ND+=1;
		}else{
			*TEMP_U+=0.5*(VX[i]*VX[i]+VY[i]*VY[i]);
			*NU+=1;
		}
	}

	*TEMP_D=*TEMP_D/(*ND);
	*TEMP_U=*TEMP_U/(*NU);
}