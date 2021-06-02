#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "MT.h"
#include <time.h>

#define Nx 30
#define Ny 15
#define PARTNUM 2
#define N (Nx*Ny)

//server
#define NSTEPS 100000000
//local test
//#define NSTEPS 200000

#define NDEVIDE 100000
#define NWRITE (NSTEPS/NDEVIDE)
#define TEMPL 4.0
#define OMEGA_0 -0.01
#define P 0.5 //透過確率
#define KIZAMI 3e-4


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


//sumnation array index

int int_array_sum(int array[N]){
//    int len;
    int i;
    int sum=0;
//    len=sizeof(array)/sizeof(int);
    for (i=0;i<N;i++){
        sum+=array[i];
    }
    return sum;
}

double double_array_sum(double array[N]){
    int i;
    double sum=0;
    for (i=0;i<N;i++){
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




int main(int argc,char *argv[]) {
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
	FILE *mini_file;
	//for animation
	FILE *rx_list;
	FILE *ry_list;
	FILE *theta_file;
    FILE *efile;
	FILE *press_file;
	FILE *temp_file;
	FILE *part_press_file;
	FILE *part_temp_file;
	FILE *all_condition_file;
	FILE *relation_heat_work;
    FILE *initarray;
    FILE *initfile;

    FILE *checkEnegy;

	//flexible dir name
	char name_kin[40];
	char name_pot[50];
	char name_tot[50];
	char name_ene[50];
	char name_dis[50];
	char name_pow[50];
	char name_ome[50];
	char name_the[50];
	char name_e[50];
	char name_mini[50];
	char name_press[50];
	char name_temp[50];
	char name_press_part[50];
	char name_temp_part[50];
	char name_all_condition[50];
	char name_relation_heat_work[50];
    char name_init[50];
    char name_initarray[50];
    char name_checkEnergy[50];
	const char* plotdata="./plotdata/modelCinit";
	const char* c_kin="/kin.dat";
	const char* c_pot="/pot.dat";
	const char* c_tot="/tot.dat";
	const char* c_ene="/energy.dat";
	const char* c_dis="/dis_p.dat";
	const char* c_pow="/pow_p.dat";
	const char* c_ome="/omega.dat";
	const char* c_the="/theta.dat";
	const char* c_e="/e_lis.dat";
	const char* c_mini="/mini.dat";
	const char* c_press="/press.dat";
	const char* c_temp="/temp.dat";
	const char* c_press_part="/press_part.dat";
	const char* c_temp_part="/temp_part.dat";
	const char* c_all_condition="/all_condition.dat";
	const char* c_relation_heat_work="/heat_and_work.dat";
    const char* c_makeinit="/init.dat";
    const char* c_initarray="/initarray.dat";
    const char* c_checkEnergy = "/checkEnergy.dat";
	//generate random seed from time
	init_genrand((unsigned)time(NULL));

	//const number
	const double rho = 0.2;
	const double dx = sqrt(1.0 / rho);
	const double dy = dx;
	const double lx = dx * Nx;
	const double ly = dy * Ny;


//コマンドライン引数から変更するためのパラメータ////
    //温度について
//	const double temp_l = TEMPL;
	const double temp_l=atof(argv[1]);

	//確率について
	double probabirity=P;
//	double probabirity=(double) 0.1*atoi(argv[1]);

	int P_dummy=10*(1-probabirity);
////////////////////////////////////////////////
	const double temp_h = 1.0;
	const double temp_p = (temp_h + temp_l) * 0.5;
	const double h = KIZAMI;
	const double h_rev=1/h;
	const double h2 = 0.5 * h * h;
	const int nout = NSTEPS / 10;
	const int t_end = NSTEPS * h;
	double w0=0;
	const int pow_num=4;

// kinds of parametor that I need
	int i, j, k, l,m,p;
	int t;
	int t_write=0;
	double ts = 0;
	double te = 0;


	//idealy thermal efficiency
	const double eta = 1.0 - (temp_h / temp_l);
	const double eta_ca = 1.0 - (sqrt(temp_h / temp_l));


    //flywheel
    const double rf=0.25*ly;//hankei
    const int rhof=pow(10,pow_num);//density of fw
    const int gamma = 0;//dumping by rotarion
    const double ia=0.5*M_PI*rhof*rf*rf*rf*rf;//inartia
	const double ria=1/ia;

    //shaft
    const double ls=pow(10,4);
    const double lsp=ls-2*rf;
 

    //displacerpiston
    int mdp = pow(10,pow_num);
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
    int mpp=pow(10,pow_num);
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

    //phai
    double dwx=0.5*lx;
	double dwx_b;
    double dwy=dpy+ls;
    double dwy_b;

    //connect pow-p to shaft
    double pwx=0.5*lx-rf;
	double pwx_b;
    double pwy=lsp+ppy0;
	double pwy_b;

    //initial arg,vel,acc
    double theta=0.0;

    double oDummy = OMEGA_0;
    double omega_0=OMEGA_0;
    double omega=omega_0;
    double alpha=0.0;
	double dtheta=0.0;
    double dpAngleTangent = cos(rf * cos(theta - M_PI / 2) / (ls + sin(theta - M_PI / 2)));
    double ppAngleTangent = cos(rf * cos(theta) / (lsp + sin(theta)));
    double dpAngleNormal = sin(rf * cos(theta - M_PI / 2) / (ls + sin(theta - M_PI / 2)));
    double ppAngleNormal = sin(rf * cos(theta) / (lsp + sin(theta)));
//    double wy_b=lscos+dpy0;

	//work and thermalefficiency
	double w=0.0;
    double wPP = 0.0;
    double wTemp = 0.0;
	double e=0.0;


	//追加
//	double one_cycle_w=0;
//	double one_cycle_q_in=0;


	//hit or through(piston)
//	int hit_piston=0;
//	int through_piston=0;


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


sprintf(name_kin,"%s%d%s",plotdata,(int) atof(argv[1]),c_kin);
sprintf(name_pot,"%s%d%s",plotdata,(int) atof(argv[1]),c_pot);
sprintf(name_tot,"%s%d%s",plotdata,(int) atof(argv[1]),c_tot);
sprintf(name_ene,"%s%d%s",plotdata,(int) atof(argv[1]),c_ene);
sprintf(name_dis,"%s%d%s",plotdata,(int) atof(argv[1]),c_dis);
sprintf(name_pow,"%s%d%s",plotdata,(int) atof(argv[1]),c_pow);
sprintf(name_ome,"%s%d%s",plotdata,(int) atof(argv[1]),c_ome);
sprintf(name_the,"%s%d%s",plotdata,(int) atof(argv[1]),c_the);
sprintf(name_e,"%s%d%s",plotdata,(int) atof(argv[1]),c_e);
sprintf(name_mini,"%s%d%s",plotdata,(int) atof(argv[1]),c_mini);
sprintf(name_press,"%s%d%s",plotdata,(int) atof(argv[1]),c_press);
sprintf(name_temp,"%s%d%s",plotdata,(int) atof(argv[1]),c_temp);
sprintf(name_press_part,"%s%d%s",plotdata,(int) atof(argv[1]),c_press_part);
sprintf(name_temp_part,"%s%d%s",plotdata,(int) atof(argv[1]),c_temp_part);
sprintf(name_all_condition,"%s%d%s",plotdata,(int) atof(argv[1]),c_all_condition);
sprintf(name_relation_heat_work,"%s%d%s",plotdata,(int) atof(argv[1]),c_relation_heat_work);
sprintf(name_init,"%s%d%s",plotdata,(int) atof(argv[1]),c_makeinit);
sprintf(name_initarray,"%s%d%s",plotdata,(int) atof(argv[1]),c_initarray);
sprintf(name_checkEnergy,"%s%d%s",plotdata,(int)atof(argv[1]),c_checkEnergy);

mini_file=fopen(name_mini,"w");
all_condition_file=fopen(name_all_condition,"w");
efile = fopen(name_e,"w");
omega_file=fopen(name_ome,"w");
theta_file=fopen(name_the,"w");
relation_heat_work=fopen(name_relation_heat_work,"w");
checkEnegy = fopen(name_checkEnergy , "w");

//double firstTotal = 0.5 * ( ia *omega_0 * omega_0 + mdp * dpv * dpv + mpp * ppv * ppv);
//-------------start mainroop------------------
	ts = omp_get_wtime();
	for (t = 1;t <= NSTEPS;t++) {

        dpy_b=dpy;	
        dwx_b=dwx;
        dwy_b=dwy;
        ppy_b=ppy;	
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
		//about disp-p
		dwx=0.5*lx+rf*sin(theta);
		dwy=ls+dpy0+rf*(1-cos(theta));
		drdx=dwx-dpx;
		drdx_b=dwx_b-dpx;
		drdy=dwy-dpy;
		drdy_b=dwy_b-dpy_b;

		//about power-p
		pwx=0.5*lx-rf*cos(theta);
		pwy=lsp+ppy0+rf*sin(theta);
		prdx=pwx-ppx;
		prdx_b=pwx_b-ppx;
		prdy=pwy-ppy;
		prdy_b=pwy_b-ppy_b;

        //サイドスラスト計算（力の接戦成分を求めている）
        dpAngleTangent = cos(rf * cos(theta - M_PI / 2) / (ls + sin(theta - M_PI / 2)));
        ppAngleTangent = cos(rf * cos(theta) / (lsp + sin(theta)));
        dpAngleNormal = sin(rf * cos(theta - M_PI / 2) / (ls + sin(theta - M_PI / 2)));
        ppAngleNormal = sin(rf * cos(theta) / (lsp + sin(theta)));

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
		dwx=0.5*lx+rf*sin(theta);
		dwy=ls+dpy0+rf*(1-cos(theta));
		pwx=0.5*lx-rf*cos(theta);
		pwy=lsp+ppy0+rf*sin(theta);
		dpv+=dpa*h*0.5+0.5*h*Ry1*rmdp;
		ppv+=ppa*h*0.5+0.5*h*Ry2*rmpp;
        oDummy += 0.5*(-Ry2*cos(theta)*h*rf)*ria;
		omega+=0.5*h*alpha+0.5*(-Ry2*cos(theta)*h*rf)*ria;
		//piston force
//		heatwall(h,N,ry,ry_b,vy,&q_in,&q_out,ppy,ppv,temp_l,temp_h,&fpp,ly,&hss,&dw,&q_debug);
//		piston_move_d(N,ry,ry_b,vy,vy_b,ay_b,h,h_rev,&q_in,&q_out,q_in_sum,q_out_sum,dpy,dpy_b,&dpv,&hit_piston,&through_piston,&fdp,temp_l,temp_h,&h1_d,momentum_d,kin_d,&down_hit,&down_through,mdp,&delta_mom_u,&delta_kin_u,probabirity,&q_debug);
//		piston_move_u(N,ry,ry_b,vy,vy_b,ay_b,h,h_rev,&q_in,&q_out,q_in_sum,q_out_sum,dpy,dpy_b,&dpv,&hit_piston,&through_piston,&fdp,temp_l,temp_h,&h1_d,momentum_u,kin_u,&up_hit,&up_through,mdp,&delta_mom_d,&delta_kin_d,probabirity,&q_debug);
		//It can be better
		ppa=fpp*rmpp;
		dpa=fdp*rmdp;
		dpv+=0.5*h*dpa;
		ppv+=0.5*h*ppa;
        if(t % 400000 == 0) {
            printf("");
        }
		alpha=-gamma*rf*omega*ria;
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
        //debug
        oDummy += 0.5*h*(-Ry2_af*cos(theta))*rf*ria;
		omega+=0.5*h*(-Ry2_af*cos(theta))*rf*ria;
        double dommga = oDummy*769657.9631818172 + 7696.774717855486;
		//boundary condition 
///////temprature and pressure onthepiston under the piston/////////
		//caliclate kin-energy
		//work and thermal efficiency
		w += gamma * rf * omega * dtheta;//仕事の合計
        wPP += (-Ry2_af*cos(theta)) * rf * dtheta;
        wTemp = 0.5 * ia * omega * omega;//運動エネルギー
		double ppEne = 0.5 * ppv * ppv * mpp;
        double dwork = dommga + w;
        double dpEnergy = 0.5 * mdp * dpv * dpv;
        double ppEnergy = 0.5 * mpp * ppv * ppv;
        double total = wTemp + dpEnergy + ppEnergy;
        if (t % nout == 0){
            printf("dpv : %lf , ppv : %lf\n" , dpv , ppv);
//            printf("firstTotal : %lf\n",firstTotal);
            printf("total : %lf\n",total);
			printf("firstEnegy : %lf  ,  now : %lf\n",0.5 * ia * omega_0 * omega_0 ,wTemp );
            printf("omega0 : %lf , omega : %lf\n", omega_0 , omega);
            printf("work : %lf , firstenergy - now : %lf\n" , w , 0.5 * ia * omega_0 * omega_0 - wTemp);
//            printf("connect PP : %lf , connect DP : %lf\n",sqrt((pwx - ppx) * (pwx - ppx) + (pwy - ppy) * (pwy - ppy)) , sqrt((dwx - dpx) * (dwx - dpx) + (dwy - dpy) * (dwy - dpy)));
//            printf("Power Piston : %lf , Displacer Piston : %lf\n\n", lsp , ls);
        }
        if (t % 10000000 == 0) {
            printf(" ");
        }
	if (t%NDEVIDE==0){
        fprintf(checkEnegy,"%lf    %lf\n",(double) t*h,total);
		fprintf(omega_file,"%lf    %lf\n",(double) t*h,w);
		fprintf(theta_file,"%lf    %lf\n",(double) t*h,0.5 * ia * omega_0 *omega_0);
	}
}
//end mainloop

te = omp_get_wtime();

//file close
fclose(all_condition_file);
fclose(mini_file);
fclose(omega_file);
fclose(theta_file);
fclose(efile);
fclose(relation_heat_work);
fclose(checkEnegy);
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
		}
	}
	for (i = 3;i < N_GY -3;i++) {
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

void force_nos(int NP, double RX[], double RY[], double AX[], double AY[], double LX, double LY, double RC2, double* POT, double* POT_IJ, double* MINI) {
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
			} else if (rxij < -0.5 * LX) {
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
			} else {
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




void heatwall(double H,int NP,double RY[],double RY_B[],double VY[],double *Q_IN,double *Q_OUT ,double PPY,double PPV,double TEMP_L,double TEMP_H,double *FPP,double LY,double *h_ss,double *d_w,double *q_debug){
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
			VY[i]=canon_b(TEMP_L);
			RY[i]+=VY[i]*(*h_ss);
			dq=0.5*(VY[i]*VY[i]-vy_l*vy_l);
            *q_debug += dq;
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
            *q_debug += dq;
			if(dq>=0.0){
				*Q_IN+=dq;
			}else{
				*Q_OUT+=dq;
			}
		}
	}
		*FPP=fpu_sum;
}




void piston_move_u(int NP,double RY[N],double RY_B[N],double VY[N],double VY_B[N],double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,double *H1_D,double Momentum_u[2],double Kin_u[2],int *k,int *j,double MDP,double *delta_mom,double *delta_kin,double PROBABIRITY,double *q_debug){
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

	for(i=0;i<N;i++){
		if(RY_B[i]>=DPY_B+0.5 && RY[i]<DPY+0.5){
            //devide step h -> h1,h2 h1
			rand_num=genrand_real3();
			d_w=RY_B[i]-(DPY_B+0.5);
			re_vel=(VY_B[i]-(*DPV))*(VY_B[i]-(*DPV));
			re_vel=sqrt(re_vel);
			h1=d_w/re_vel;
			h_ss=fabs(H-h1);
			RY[i]=RY_B[i]+VY_B[i]*h1;
            DPY=DPY_B+(*DPV)*h1;
            vy_l=VY[i];
            if (rand_num < (1 - PROBABIRITY)) {         
                *j+=1;
                *THROUGH_PISTON+=1;
                VY[i]=-canon_b(TEMP_L);
                dd=VY[i]*h_ss;
                RY[i]=RY[i]+dd;
                fdu=-(VY[i]-vy_l)*H_REV;
                fdu_sum+=fdu;
                dq=0.5*(VY[i]*VY[i]-vy_l*vy_l);
                *q_debug += dq;
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
	
void piston_move_d(int NP,double RY[N],double RY_B[N],double VY[N],double VY_B[N],double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,double *H1_D,double Momentum_d[2],double Kin_d[2],int *kk,int *jj,double MDP,double *delta_mom,double *delta_kin,double PROBABIRITY,double *q_debug){
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
    for(i=0;i<N;i++){
    	if(RY_B[i]<=DPY_B-0.5 && RY[i]>DPY-0.5){
			z=genrand_real3();
			d_w=(DPY_B-0.5)-RY_B[i];
			re_vel=(VY_B[i]-(*DPV))*(VY_B[i]-(*DPV));
			re_vel=sqrt(re_vel);
			h1=d_w/re_vel;
			h_ss=fabs(H-h1);
			RY[i]=RY_B[i]+VY_B[i]*h1;
			DPY=DPY_B+(*DPV)*h1;
            vy_l=VY[i];         
			if (z<(1-PROBABIRITY))
                *jj+=1;
				*THROUGH_PISTON+=1;
                VY[i]=canon_b(TEMP_H);
				dd=VY[i]*h_ss;
				RY[i]=RY[i]+dd;
				fdu=-(VY[i]-vy_l)*H_REV;
				fdu_sum+=fdu;
				dq=0.5*(VY[i] * VY[i] - vy_l * vy_l);
                *q_debug += dq;
                if(dq>=0.0){
                    *Q_IN+=dq;
                }else{
                    *Q_OUT+=dq;
                }                     
        }
    }
	*FDP=fdd_sum+fdu_sum;
}


void velocity(int NP,double RY[],double VX[],double VY[],double DPY,double *TEMP_D,double *TEMP_U,int *ND,int*NU ,int N_U_list[N],double *TEMP){
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
	*TEMP=*TEMP/N;
	*TEMP_D=*TEMP_D/(*ND);
	*TEMP_U=*TEMP_U/(*NU);
}


void partition(double RY[],double VX[] , double VY[] ,double PPY,double TEMP_PART_LIST[NWRITE][PARTNUM],double PRESS_PART_LIST[NWRITE][PARTNUM],int T_WRITE){
	int i,j,k;
	double temp;
	double re_partnum=(double) 1/PARTNUM;
	double part,part_plus;
	double pow_ppy=(PPY*re_partnum)*(PPY*re_partnum);
	for(i=0;i<PARTNUM;i++){
		k=0;
		part=((double) (i/PARTNUM))*PPY;
		part_plus=((double) (i+1)/PARTNUM)*PPY;
		for(j=0;j<N;j++){
			if(part<=RY[j] && part_plus>RY[j]){
				k+=1;
				TEMP_PART_LIST[T_WRITE-1][i]+=0.5*(VX[j]*VX[j]+VY[j]*VY[j]);
				}
			}
		TEMP_PART_LIST[T_WRITE-1][i]=TEMP_PART_LIST[T_WRITE-1][i]/k;
		temp=TEMP_PART_LIST[T_WRITE-1][i];
		PRESS_PART_LIST[T_WRITE-1][i]=(k/pow_ppy)*temp;
	
	}
}
