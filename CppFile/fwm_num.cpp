#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "MT.h"
#include <time.h>



#define Nx 10
#define Ny 20
#define PARTNUM 2
#define NUM (Nx*Ny)

/*本番環境*/
//#define NSTEPS 1500000000
/*ローカルサーバ*/
#define NSTEPS 200000

#define NDEVIDE 100000
#define NWRITE (NSTEPS/NDEVIDE)
#define TEMPL 4.0
#define OMEGA_0 -0.01
#define P 0.5
#define KIZAMI 1e-3

/*関数プロトタイプ宣言*/
void gmap_create(int NP, double RX[], double RY[], double L_GX, double L_GY, int N_GX, int N_GY, int NEIGHBOR_ROW[], int NEIGHBOR_COL[], int NEITGHBOR_LEN, int PAIRLIST[][10], double LX, int **G_MAP);

void force(int NP, double RX[], double RY[], double AX[], double AY[], double LX, double LY, int PAIRLIST[][10], double RC2, double* POT, double* POT_IJ, double* MINI);

void heatwall(double H,int NP,double RY[],double RY_B[],double VY[],double *Q_IN,double *Q_OUT ,double PPY,double PPV,double TEMP_L,double TEMP_H,double *FPP,double LY,double *h_ss,double *d_w);

void piston_move_u(int NP,double RY[],double RY_B[],double VY[],double VY_B[],double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,double *H1_D,double Momentum_u[2],double Kin_u[2],int *k,int *j,double MDP,double *delta_mom,double *delta_kin,double PROBABIRITY);

void piston_move_d(int NP,double RY[],double RY_B[],double VY[],double VY_B[],double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,double *H1_D,double Momentum_d[2],double Kin_d[2],int *kk,int *jj,double MDP,double *delta_mom,double *delta_kin,double PROBABIRITY);

void velocity(int NP,double RY[],double VX[],double VY[],double DPY,double *TEMP_D,double *TEMP_U,int *ND,int*NU ,int N_U_list[],double *TEMP);

void partition(double RY[],double VX[] , double VY[] ,double PPY,double TEMP_PART_LIST[NWRITE][PARTNUM],double PRESS_PART_LIST[NWRITE][PARTNUM],int T_WRITE);


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


/*配列要素の合計を出す関数*/
template <typename T> T array_sum(T array[]){
	T sum=0;
	for(int i=0;i<NUM;i++){
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


/*ボルツマン分布生成*/
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


/*境界条件*/
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

	char plotdata[20]="./plotdata/m2a_num";

	char c_kin[20]="/kin.dat";
	char c_pot[20]="/pot.dat";
	char c_tot[20]="/tot.dat";
	char c_ene[20]="/energy.dat";
	char c_dis[20]="/dis_p.dat";
	char c_pow[20]="/pow_p.dat";
	char c_ome[20]="/omega.dat";
	char c_the[20]="/theta.dat";
	char c_e[20]="/e_lis.dat";
	char c_mini[20]="/mini.dat";
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


/*コマンドライン引数から変更するためのパラメータ*/
    /*温度について*/
	const double temp_l = TEMPL;
//	const double temp_l=atof(argv[1]);

/*確率について*/
	double probabirity=P;
//	double probabirity=(double) 0.1*atoi(argv[1]);

	int P_dummy=10*(1-probabirity);

/*粒子数を変数に保存*/
    const int N=NUM;

/*上壁の温度*/
	const double temp_h = 1.0;

/*刻み幅とその逆数,乗数*/
	const double h = KIZAMI;
	const double h_rev=1/h;
	const double h2 = 0.5 * h * h;

/*出力頻度と終了時刻*/
	const int nout = NSTEPS / 10;
	const int t_end = NSTEPS * h;

/*運動エネルギー、ポテンシャル、全エネルギー*/
	double kin0=0;
	double pot = 0;
	double pot_ij = 0;
	double total_kin=0;
	double total_e=0;
	const int POW_NUM=4;

/*ディスプレーサに当たった粒子、すり抜けた粒子*/
	int up_hit=0;
	int down_hit=0;
    int up_through=0;
    int down_through=0;
	int hit_piston=0;
	int through_piston=0;

/*モーメントとエネルギーのチェック（ディスプレーサ）*/
	double momentum_u[2]={0.0,0.0};
	double kin_u[2]={0.0,0.0};
	double momentum_d[2]={0.0,0.0};
	double kin_d[2]={0.0,0.0};
	double delta_mom_u=0.0;
	double delta_mom_d=0.0;
	double delta_kin_u=0.0;
	double delta_kin_d=0.0;

/*時刻に対する各パラメータ*/
	int i, j, k, l,m,p;
	int t;
	int t_write=0;
	double ts = 0;
	double te = 0;

/*粒子の位置、速度、加速度*/
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
	double sum_vx=0.0;
	double sum_vy=0.0;

/*グリッド作成パラメータ*/
	int bo=(lx/0.4);
	double l_gx_d=lx/bo;
	const double l_gx = l_gx_d;
	const double l_gy = l_gx;
	const int n_gx = ceil(lx / l_gx) +6 ;//袖領域のため
	const int n_gy = ceil((ly / l_gy)+(0.25*lx)/l_gy)+6;//袖領域のため
	const int n_g_all = n_gx * n_gy;
/*
	const double l_gx = 0.448;
	const double l_gy = l_gx;
	const int n_gx = ceil(lx / l_gx) + 6;//袖領域のため
	const int n_gy = ceil(ly / l_gy) + 6 + ceil(0.25*lx/l_gx);//袖領域のため
	const int n_g_all = n_gx * n_gy;
*/

/*グリッド配列作成*/
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


							//									/////////////////////	
						  /////									//	               //
						/////////								//	   -3-3-3      //
							//									//	   -2-2-2      //
							//									//	   -1-1-1      //
							//									//	    0 0 0      //
							////////////////////////////////	//      1 1 1 1    //
																//	    2 2 2 2    //
																//	    3 3 3 3    //
																//                 // 
																/////////////////////	

/*ペアリスト作成*/
	int pairlist[N][10];
	int pairlist_innit[N][10];

/*理想的な熱効率*/
	const double eta = 1.0 - (temp_h / temp_l);
	const double eta_ca = 1.0 - (sqrt(temp_h / temp_l));

/*力の計算時に使うパラメータ*/
	const double rc = pow(2.0, (double)1 / 6);
	double rc2 = rc * rc;
	double mini = 1000.0;
	double mini_b=1000.0;

/*フライホイール*/
    const double rf=10;//hankei
    const int rhof=2*pow(10,POW_NUM);//density of fw
    const int gamma=100;//dumping by rotarion
//    const double mf=M_PI*rf*rf*rhof;
    const double ia=0.5*M_PI*rhof*rf*rf*rf*rf;//inartia
	const double ria=1/ia;

/*シャフト*/
    const double ls=pow(10,4);
    const double lsp=ls-2*rf;

/*熱壁*/
    double q_in=0.0;
    double q_out=0.0;
    double q_in_sum=0.0;
    double q_out_sum=0.0;

/*ディスプレーサピストン*/
    int mdp=2*pow(10,POW_NUM);
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
    const double mu1=ia*mdp/(ia+mdp);//mass of calculate

/*パワーピストン*/
    int mpp=2*pow(10,POW_NUM);
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
    const double mu2=ia*mpp/(ia+mpp);



/*ディスプレーサピストン-シャフト*/
    double lscos=ls*(sqrt(ls*ls-rf*rf)/ls);
    //phai
    double dwx=0.5*lx;
	double dwx_b;
    double dwy=dpy+ls;
    double dwy_b;

/*パワーピストン-シャフト*/
    double pwx=0.5*lx-rf;
	double pwx_b;
    double pwy=lsp+ppy0;
	double pwy_b;

/*回転角、初期角速度、角加速度*/
    double theta=0.0;
    double omega_0=OMEGA_0;
    double omega=omega_0;
    double alpha=0.0;
	double dtheta=0.0;
    double wy_b=lscos+dpy0;

/*仕事、熱効率*/
	double w=0.0;
	double e=0.0;

/*系の上下別温度、圧力*/
	double temp_d=0.0;
	double temp_u=0.0;
	double temp_all=0.0;
	double temp_lis[NWRITE];
	double press_d=0.0;
	double press_u=0.0;
	double press_all=0.0;
	double press_lis[NWRITE];

/*系の上下別粒子数*/
	int N_D=0;
	int N_U=0;
	int N_D_b=0;
	int N_U_b=0;
    int N_U_init[N];
    int N_U_l[N];
	int N_U_l_b_lis[N];
    int N_U_l_b;
    int N_U_init_sum=0;
    int N_U_l_sum=0;

/*拘束パラメータ*/
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

/*パワーピストンにはたらく力*/
	double fp_sum=0.0;
	double fd_sum=0.0;

/*ディスプレーサピストンにはたらく力*/
	double fpp=0.0;
	double fdp=0.0;

/*デバッグに使うパラメータ*/
	double hss=0.0;
	double dw=0.0;
    double h1_d=0.0;

/*粒子の位置情報読み込み*/
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

/*大きい配列*/
	double t_lis[NWRITE];

/*回転角ごとに出力するためのパラメータ*/
	const double ack_rot=0.125*M_PI;
    const double ack_circle =2*M_PI;
	double rot=0.0;
    double circle=0.0;
	int rot_num=0;
    int b_num=0;
    int circle_num=0;

/*系内を区切った時の部分区間の温度、圧力*/
	double temp_part_lis[NWRITE][PARTNUM];
	double press_part_lis[NWRITE][PARTNUM];

/*部分圧力、温度の初期化*/
	for(i<0;i<NWRITE;i++){
		for(j=0;j<PARTNUM;j++){
			temp_part_lis[i][j]=0.0;
			press_part_lis[i][j]=0.0;
		}
	}

/*グリッドマップ、ペアリスト初期化*/	
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

/*粒子の初期速度*/
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

/*粒子の初期位置、初期各加速度*/
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

/*粒子数リスト初期化*/
    for(i=0;i<N;i++){
        N_U_init[i]=-1;
        N_U_l[i]=-1;
    }

/*初期運動エネルギー*/
	for (i = 0;i < N;i++) {
		kin0 = kin0 + (vx[i] * vx[i] + vy[i] * vy[i]);
	}
	kin0 = kin0 * 0.5;

/*初期ペアリスト作成*/
	const int neighbor_len = sizeof(neighbor_list_row) / sizeof(int);
	gmap_create(N, rx, ry, l_gx, l_gy, n_gx, n_gy, neighbor_list_row, neighbor_list_col, neighbor_len, pairlist, lx, g_map);

/*ペアリスト作成確認*/
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
	printf("kin0:%lf\n",kin0);

velocity(N,ry,vx,vy,dpy,&temp_d,&temp_u,&N_D_b,&N_U_b,N_U_init,&temp_all);
velocity(N,ry,vx,vy,dpy,&temp_d,&temp_u,&N_D,&N_U,N_U_l,&temp_all);
N_U_init_sum=array_sum(N_U_init);
//rx_list=fopen("./plotdata/rx_list.dat","w");
//ry_list=fopen("./plotdata/ry_list.dat","w");
//dp=fopen("./plotdata/pack2_wt/dp.dat","w");
//py=fopen("./plotdata/pack2_wt/pp.dat","w");

/*ファイルパス作成*/
sprintf(name_kin,"%s%d%s",plotdata,(int) N,c_kin);
sprintf(name_pot,"%s%d%s",plotdata,(int) N,c_pot);
sprintf(name_tot,"%s%d%s",plotdata,(int) N,c_tot);
sprintf(name_ene,"%s%d%s",plotdata,(int) N,c_ene);
sprintf(name_dis,"%s%d%s",plotdata,(int) N,c_dis);
sprintf(name_pow,"%s%d%s",plotdata,(int) N,c_pow);
sprintf(name_ome,"%s%d%s",plotdata,(int) N,c_ome);
sprintf(name_the,"%s%d%s",plotdata,(int) N,c_the);
sprintf(name_e,"%s%d%s",plotdata,(int) N,c_e);
sprintf(name_mini,"%s%d%s",plotdata,(int) N,c_mini);
sprintf(name_press,"%s%d%s",plotdata,(int) N,c_press);
sprintf(name_temp,"%s%d%s",plotdata,(int) N,c_temp);
sprintf(name_press_part,"%s%d%s",plotdata,(int) N,c_press_part);
sprintf(name_temp_part,"%s%d%s",plotdata,(int) N,c_temp_part);
sprintf(name_all_condition,"%s%d%s",plotdata,(int) N,c_all_condition);
mini_file=fopen(name_mini,"w");
all_condition_file=fopen(name_all_condition,"w");
efile = fopen(name_e,"w");
omega_file=fopen(name_ome,"w");
theta_file=fopen(name_the,"w");

/*メインループはじまり*/
	ts = omp_get_wtime();
	for (t = 1;t <= NSTEPS;t++)
	{
		/*一つ前のステップの保存*/
		for(i=0;i<N;i++){
			rx_b[i]=rx[i];
			ry_b[i]=ry[i];
			vy_b[i]=vy[i];
			ay_b[i]=ay[i];
			N_U_l_b_lis[i]=N_U_l[i];
		}
			dpy_b=dpy;	
			dwx_b=dwx;
			dwy_b=dwy;
			ppy_b=ppy;	
			N_U_l_b=N_U_l_sum;

/*
		if(t%10==0){
			rx_list=fopen("./plotdata/rx_list.dat","w");
			ry_list=fopen("./plotdata/ry_list.dat","w");
			for (i=0;i<N;i++){
				fprintf(rx_list,"%lf    ",rx[i]);
				fprintf(ry_list,"%lf    ",ry[i]);
			}
			fclose(rx_list);
			fclose(ry_list);
		}
*/

/*粒子最短距離更新*/
		if (mini_b>mini){
			fprintf(mini_file,"%lf\n",mini);
			mini_b=mini;
		}



/*初期化*/
		total_kin = 0.0;
		total_e = 0.0;
		if ((t-1)%NDEVIDE == 0){
			t_write+=1;
			t_lis[t_write-1] = t * h;
			//圧力などの定義
			partition(ry,vx,vy,ppy,temp_part_lis,press_part_lis,t_write);
		}



/*ベルレ法一回目*/
		for (i = 0;i < N;i++) {
			rx[i] = rx[i] + vx[i] * h + ax[i] * h2;
			ry[i] = ry[i] + vy[i] * h + ay[i] * h2;
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}

/*-1クリア*/
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

/*グリッドマップ作成*/
		gmap_create(N, rx, ry, l_gx, l_gy, n_gx, n_gy, neighbor_list_row, neighbor_list_col, neighbor_len, pairlist, lx, g_map);

/*力の計算*/		
		force(N, rx, ry, ax, ay, lx, ly, pairlist, rc2, &pot, &pot_ij, &mini);

/*ベルレ法二回目*/
		for (i = 0;i < N;i++) {
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}

/*速度の合計を計算*/
		sum_vx=array_sum(vx);
		sum_vy=array_sum(vy);

/*ピストンの速度ベルレ法*/
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

/*シャフトーピストンの運動方程式*/
/*拘束力学系計算(ディスプレーサーピストン)*/
		dwx=0.5*lx+rf*sin(theta);
		dwy=ls+dpy0+rf*(1-cos(theta));
		drdx=dwx-dpx;
		drdx_b=dwx_b-dpx;
		drdy=dwy-dpy;
		drdy_b=dwy_b-dpy_b;

/*拘束力学系計算(パワーピストン)*/
		pwx=0.5*lx-rf*cos(theta);
		pwy=lsp+ppy0+rf*sin(theta);
		prdx=pwx-ppx;
		prdx_b=pwx_b-ppx;
		prdy=pwy-ppy;
		prdy_b=pwy_b-ppy_b;

/*拘束力計算、回転角更新*/
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

/*ベルレ法更新*/
		dwx=0.5*lx+rf*sin(theta);
		dwy=ls+dpy0+rf*(1-cos(theta));
		pwx=0.5*lx-rf*cos(theta);
		pwy=lsp+ppy0+rf*sin(theta);
		dpv+=dpa*h*0.5+0.5*h*Ry1*rmdp;
		ppv+=ppa*h*0.5+0.5*h*Ry2*rmpp;
		omega+=0.5*h*alpha+0.5*(-Ry2*cos(theta)*h*rf)*ria;

/*熱壁、パワーピストン、ディスプレーさーピストン計算*/
		heatwall(h,N,ry,ry_b,vy,&q_in,&q_out,ppy,ppv,temp_l,temp_h,&fpp,ly,&hss,&dw);
		piston_move_d(N,ry,ry_b,vy,vy_b,ay_b,h,h_rev,&q_in,&q_out,q_in_sum,q_out_sum,dpy,dpy_b,&dpv,&hit_piston,&through_piston,&fdp,temp_l,temp_h,&h1_d,momentum_d,kin_d,&down_hit,&down_through,mdp,&delta_mom_u,&delta_kin_u,probabirity);
		piston_move_u(N,ry,ry_b,vy,vy_b,ay_b,h,h_rev,&q_in,&q_out,q_in_sum,q_out_sum,dpy,dpy_b,&dpv,&hit_piston,&through_piston,&fdp,temp_l,temp_h,&h1_d,momentum_u,kin_u,&up_hit,&up_through,mdp,&delta_mom_d,&delta_kin_d,probabirity);

/*拘束力を含めたベルレ法の更新*/
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

/*周期境界条件適用*/
		boundary(N, rx, lx);

/*温度圧力計算*/
		velocity(N,ry,vx,vy,dpy,&temp_d,&temp_u,&N_D,&N_U,N_U_l,&temp_all);

		press_all=temp_all*N/(lx*ppy);
		press_d=temp_d*N_D/(lx*dpy);
		press_u=temp_u*N_U/(lx*(ppy-dpy));
        N_U_l_sum=array_sum(N_U_l);


/*運動エネルギー*/
		for (i = 0;i < N;i++) {
			total_kin += (vx[i] * vx[i] + vy[i] * vy[i]);
		}
		total_kin = total_kin * 0.5;

/*全エネルギー*/
		total_e = total_kin + pot;

/*計算状況出力*/
			if (t % nout == 0){
				printf("Time:%lf,Kinetic Energy:%lf,Potential Energy:%lf,Total Energy:%lf,momentum_vx:%lf,momentum_vy:%lf\n", t * h, total_kin, pot, total_e,sum_vx,sum_vy);
				}

/*仕事と熱効率計算*/
		w+=gamma*rf*omega*dtheta;
		e=w/(0.5*ia*omega_0*omega_0+q_in);

/*PV線図に必要な回転状況計算*/
		rot=fabs(theta)-rot_num*ack_rot;
		if (rot>ack_rot){
			rot-=fabs(ack_rot);
			fprintf(all_condition_file,"%lf    %lf    %lf    %lf    %lf    %lf\n",t*h,press_all,temp_all,ppy*lx,ppy,dpy);
			rot_num+=1;
		}

/*回転数計算*/
        circle=fabs(theta)-rot_num*ack_circle;
		if (circle>ack_circle){
			circle-=fabs(ack_circle);
			circle_num+=1;
		}

/*ファイル出力*/
	if (t%NDEVIDE==0){
		fprintf(omega_file,"%lf    %lf\n",(double) t*h,omega);
		fprintf(theta_file,"%lf    %lf\n",(double) t*h,theta);
		fprintf(efile,"%lf    %lf\n",(double) t*h,e);
	}

/*
		//into the big array
		if ((t-1)%NDEVIDE == 0)
		{
			e_lis[t_write-1]=e;
			omega_lis[t_write-1]=omega;
			theta_lis[t_write-1]=theta;
			press_lis[t_write-1]=press_all;
			temp_lis[t_write-1]=temp_all;


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
	
		}
*/	

}

/*メインループおわり*/
te = omp_get_wtime();
printf("time cost is %lf seconds\n", te - ts);
printf("process time is:%lf\n",te-ts);
printf("hit is :%d,through is %d\n",hit_piston,through_piston);
printf("up_hit:%d,down_hit:%d\n",up_hit,down_hit);
printf("up_through:%d,down_through:%d\n",up_through,down_through);
printf("rotation num is : %d",circle_num);

/*メモリ開放*/
for(i=0;i<n_gy;i++){
	free(g_map[i]);
}
free(g_map);

/*ファイル出力*/
//press_file=fopen(name_press,"w");
//temp_file=fopen(name_temp,"w");
//part_press_file=fopen(name_press_part,"w");
//part_temp_file=fopen(name_temp_part,"w");
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

/*ファイル閉じる*/
fclose(all_condition_file);
fclose(mini_file);
fclose(omega_file);
fclose(theta_file);
fclose(efile);
//fclose(press_file);
//fclose(temp_file);
//fclose(part_press_file);
//fclose(part_temp_file);
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
//debugfile


//lxly=fopen("./plotdata/pack2_wt/lxly.dat","w");
//fprintf(lxly,"%lf    %lf",lx,ppy_max);
//fclose(lxly);




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




void piston_move_u(int NP,double RY[],double RY_B[],double VY[],double VY_B[],double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,double *H1_D,double Momentum_u[2],double Kin_u[2],int *k,int *j,double MDP,double *delta_mom,double *delta_kin,double PROBABIRITY){
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
           if (rand_num<(1-PROBABIRITY)){
				*HIT_PISTON+=1;
    			*k+=1;
				Momentum_u[0]=VY[i]+MDP*(*DPV);
				Kin_u[0]=0.5*VY[i]*VY[i]+0.5*MDP*(*DPV)*(*DPV);

				VY[i]=-((MDP-1)/(MDP+1))*VY_B[i]+((2*MDP)/(MDP+1))*(*DPV);
              *DPV=((MDP-1)/(MDP+1))*(*DPV)+(2/(MDP+1))*vy_l;

//				VY[i]=-VY_B[i]+2*(*DPV);

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
                VY[i]=-canon_b(TEMP_L);
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
	
void piston_move_d(int NP,double RY[],double RY_B[],double VY[],double VY_B[],double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,double *H1_D,double Momentum_d[2],double Kin_d[2],int *kk,int *jj,double MDP,double *delta_mom,double *delta_kin,double PROBABIRITY){
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
    for(i=0;i<NUM;i++){
    	if(RY_B[i]<=DPY_B-0.5 && RY[i]>DPY-0.5){

			z=genrand_real3();
			d_w=(DPY_B-0.5)-RY_B[i];
			re_vel=(VY_B[i]-(*DPV))*(VY_B[i]-(*DPV));
			re_vel=sqrt(re_vel);
			h1=d_w/re_vel;
			h_ss=fabs(H-h1);

			RY[i]=RY_B[i]+VY_B[i]*h1;
			DPY=DPY_B+(*DPV)*h1;
			if (z<(1-PROBABIRITY)){
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
	*TEMP=*TEMP/NUM;
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
		for(j=0;j<NUM;j++){
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
