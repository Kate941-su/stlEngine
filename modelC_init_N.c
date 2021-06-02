/////////////////////////////////////////////////////////
//コマンドライン引数はなし、直接ヘッダから粒子数をいじる//
/////////////////////////////////////////////////////////
#include "modelC_init.h"
#include <omp.h>
#include <time.h>
#define TEMPL 4.0
#define OMEGA_0 -0.01
#define P 0.5 //透過確率
#define KIZAMI 3e-4

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
    FILE *metadata;

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
	const char* plotdata="./plotdata/modelCinitN";
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
	const double temp_l=TEMPL;
	//確率について
	double probabirity=P;
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
	double kin0=0;
	double pot = 0;
	double pot_ij = 0;
	double total_kin=0;
	double total_e=0;
	const int pow_num=4;
	int up_hit=0;
	int down_hit=0;
    int up_through=0;
    int down_through=0;
// kinds of parametor that I need
	int i, j, k, l,m,p;
	int t;
	int t_write=0;
	double ts = 0;
	double te = 0;
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
	int pairlist[N][10];
	int pairlist_innit[N][10];

	//排他的グリッド高速化法　
	int bo=(lx/0.4);
	double l_gx_d=lx/bo;
	const double l_gx = l_gx_d;
	const double l_gy = l_gx;
	const int n_gx = ceil(lx / l_gx) +6 ;//袖領域のため
	const int n_gy = ceil((ly / l_gy)+(0.25*lx)/l_gy)+6;//袖領域のため
	const int n_g_all = n_gx * n_gy;
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
	double mini_b=1000.0;
	//debug parameter
	double FX=0.0;
	double FY=0.0;

    //flywheel
    const double rf=0.25*ly;//hankei
    const int rhof=pow(10,pow_num);//density of fw
    const int gamma=100;//dumping by rotarion
    const double ia=0.5*M_PI*rhof*rf*rf*rf*rf;//inartia
	const double ria=1/ia;

    //shaft
    const double ls=pow(10,4);
    const double lsp=ls-2*rf;

	//heatwall
    double q_debug = 0.0;
    double q_in=0.0;
    double q_out=0.0;
    double q_in_sum=0.0;
    double q_out_sum=0.0;

    //ディスプレーサピストンのパラメータ
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

    //パワーピストンのパラメータ
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

    //シャフト
    double lscos=ls*(sqrt(ls*ls-rf*rf)/ls);

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
    double omega_0=OMEGA_0;
    double omega=omega_0;
    double alpha=0.0;
	double dtheta=0.0;
    double wy_b=lscos+dpy0;

	//work and thermalefficiency
	double w=0.0;
    double wTemp = 0.0;
	double e=0.0;


	//追加
	double one_cycle_w=0;
	double one_cycle_q_in=0;


	//hit or through(piston)
	int hit_piston=0;
	int through_piston=0;

	//tempreture and pressure
	double temp_d=0.0;
	double temp_u=0.0;
	double temp_all=0.0;
	double press_d=0.0;
	double press_u=0.0;
	double press_all=0.0;
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

	//11/30追加コンテンツ
	const double ack_rot=0.125*M_PI;
    const double ack_circle =2*M_PI;
	double rot=0.0;
    double circle=0.0;
	int rot_num=0;
    int b_num=0;
    int circle_num=0;

	double temp_part_lis[NWRITE][PARTNUM];
	double press_part_lis[NWRITE][PARTNUM];

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
    //initialize N_list
    for(i=0;i<N;i++){
        N_U_init[i]=-1;
        N_U_l[i]=-1;
    }
	for (i = 0;i < N;i++) {
		kin0 = kin0 + (vx[i] * vx[i] + vy[i] * vy[i]);
	}
	kin0 = kin0 * 0.5;
    w0 = 0.5 * ia * omega_0 * omega_0;
    double dispEne = 0.5 * mdp * dpv *dpv;//ディスプレーサピストンの運動エネルギー
    double powerEne = 0.5 * mpp * ppv * ppv;//パワーピストンの運動エネルギー
    double dispEne0 = dispEne;//ディスプレーサピストンの運動エネルギー
    double powerEne0 = powerEne;//パワーピストンの運動エネルギー
    ppy_b=ppy0;
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
		for(j=0;j<PARTNUM;j++){
			temp_part_lis[i][j]=0.0;
			press_part_lis[i][j]=0.0;
		}
	}

    printf("particlenum:%d\n",l);
    printf("kin0:%lf\n",kin0);
    calculateTemleture(N,ry,vx,vy,dpy,&temp_d,&temp_u,&N_D_b,&N_U_b,N_U_init,&temp_all);
    calculateTemleture(N,ry,vx,vy,dpy,&temp_d,&temp_u,&N_D,&N_U,N_U_l,&temp_all);
    N_U_init_sum=int_array_sum(N_U_init);

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
    sprintf(name_relation_heat_work,"%s%d%s",plotdata,(int) N,c_relation_heat_work);
    sprintf(name_init,"%s%d%s",plotdata,(int) N,c_makeinit);
    sprintf(name_initarray,"%s%d%s",plotdata,(int) N,c_initarray);

    mini_file=fopen(name_mini,"w");
    all_condition_file=fopen(name_all_condition,"w");
    efile = fopen(name_e,"w");
    omega_file=fopen(name_ome,"w");
    theta_file=fopen(name_the,"w");
    relation_heat_work=fopen(name_relation_heat_work,"w");
//-------------start mainroop------------------
	ts = omp_get_wtime();
	for (t = 1;t <= NSTEPS;t++) {
		//save before data
		for(i=0;i<N;i++){
			rx_b[i]=rx[i];
			ry_b[i]=ry[i];
			vy_b[i]=vy[i];
			ay_b[i]=ay[i];
			N_U_l_b_lis[i]=N_U_l[i];
		}
		if (mini_b>mini){
			fprintf(mini_file,"%lf\n",mini);
			mini_b=mini;
		}
        dpy_b=dpy;	
        dwx_b=dwx;
        dwy_b=dwy;
        ppy_b=ppy;	
        N_U_l_b=N_U_l_sum;
        double omega_b = omega;
		//initialize
		total_kin = 0.0;
		total_e = 0.0;
		if ((t-1)%NDEVIDE == 0){
			t_write+=1;
			//圧力などの定義
			calcTempAndPress(ry,vx,vy,ppy,temp_part_lis,press_part_lis,t_write);
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
		for (i = 0; i < N; i++) {
			for (j = 0; j < 10; j++) {
				pairlist[i][j] = -1;
			}
		}
		gmap_create(N, rx, ry, l_gx, l_gy, n_gx, n_gy, neighbor_list_row, neighbor_list_col, neighbor_len, pairlist, lx, g_map);
		//cliculate force		
		force(N, rx, ry, ax, ay, lx, ly, pairlist, rc2, &pot, &pot_ij, &mini);
		//second verlet
		for (i = 0;i < N;i++) {
			vx[i] = vx[i] + ax[i] * h * 0.5;
			vy[i] = vy[i] + ay[i] * h * 0.5;
		}

		sum_vx=double_array_sum(vx);
		sum_vy=double_array_sum(vy);

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
		omega+=0.5*h*alpha+0.5*(-Ry2*cos(theta)*h*rf)*ria;
		//piston force
		heatwall(h,N,ry,ry_b,vy,&q_in,&q_out,ppy,ppv,temp_l,temp_h,&fpp,ly,&hss,&dw,&q_debug);
		piston_move_d(N,ry,ry_b,vy,vy_b,ay_b,h,h_rev,&q_in,&q_out,q_in_sum,q_out_sum,
                    dpy,dpy_b,&dpv,&hit_piston,&through_piston,&fdp,temp_l,temp_h,
                    &h1_d,&down_hit,&down_through,mdp,probabirity,&q_debug);
		piston_move_u(N,ry,ry_b,vy,vy_b,ay_b,h,h_rev,&q_in,&q_out,q_in_sum,q_out_sum,
                    dpy,dpy_b,&dpv,&hit_piston,&through_piston,&fdp,temp_l,temp_h,
                    &h1_d,&up_hit,&up_through,mdp,probabirity,&q_debug);
		//It can be better
		ppa=fpp*rmpp;
		dpa=fdp*rmdp;
		dpv+=0.5*h*dpa;
		ppv+=0.5*h*ppa;
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
		omega+=0.5*h*(-Ry2_af*cos(theta))*rf*ria;
		//boundary condition 
		boundary(N, rx, lx);
///////temprature and pressure onthepiston under the piston/////////
		calculateTemleture(N,ry,vx,vy,dpy,&temp_d,&temp_u,&N_D,&N_U,N_U_l,&temp_all);
		press_all=temp_all*N/(lx*ppy);
		press_d=temp_d*N_D/(lx*dpy);
		press_u=temp_u*N_U/(lx*(ppy-dpy));
        N_U_l_sum=int_array_sum(N_U_l);
		//caliclate kin-energy
		for (i = 0;i < N;i++) {
			total_kin += (vx[i] * vx[i] + vy[i] * vy[i]);
		}
		total_kin = total_kin * 0.5;
		total_e = total_kin + pot;

		if (t % nout == 0){
			printf("Time:%lf,Kinetic Energy:%lf,Potential Energy:%lf,Total Energy:%lf,Qin:%lf,Qout:%lf,firstEne:%lf,difference:%lf\n", 
            t * h, total_kin, pot, total_e,q_in,q_out,kin0,kin0 + (q_in + q_out) - (total_e));
            printf("calcwork: %lf\n",w);
        }
		//work and thermal efficiency
		w += gamma * rf * omega * dtheta;//仕事の合計
        wTemp = 0.5 * ia * omega * omega;//運動エネルギー
        dispEne = 0.5 * mdp * dpv * dpv;//DPのエネルギー
        powerEne = 0.5 * mpp * ppv * ppv;//PPのエネルギー

		rot=fabs(theta)-rot_num*ack_rot;
        circle=fabs(theta)-circle_num*ack_circle;
		if (rot>ack_rot){
			rot-=fabs(ack_rot);
			fprintf(all_condition_file,"%lf    %lf    %lf    %lf    %lf    %lf\n",t*h,press_all,temp_all,ppy*lx,ppy,dpy);
			rot_num+=1;
        }
		if (circle>ack_circle){
			circle-=fabs(ack_circle);
			circle_num+=1;
			one_cycle_q_in=q_in-one_cycle_q_in;
			one_cycle_w=w-one_cycle_w;
			e=one_cycle_w/one_cycle_q_in;
			fprintf(efile,"%lf    %lf\n",(double) t*h,e);
		}

	if (t%NDEVIDE==0){
		fprintf(omega_file,"%lf    %lf\n",(double) t*h,omega);
		fprintf(theta_file,"%lf    %lf\n",(double) t*h,theta);
		fprintf(relation_heat_work,"%lf    %lf    %lf    %lf    %lf\n",((double) t*h),q_in,q_out,(q_in+q_out),w);
	}
}
//end mainloop

te = omp_get_wtime();
e=w/q_in;//熱効率（仮）
printf("time cost is %lf seconds\n", te - ts);
printf("hit is :%d,through is %d\n",hit_piston,through_piston);
printf("up_hit:%d,down_hit:%d\n",up_hit,down_hit);
printf("up_through:%d,down_through:%d\n",up_through,down_through);
printf("rotation num is : %d\n",circle_num);
//free gmap
for(i=0;i<n_gy;i++){
	free(g_map[i]);
}
free(g_map);


double makeinit_array[]={dpx,dpy,ppx,ppy,dpy_b,ppy_b,dpv,ppv,ppa,dpa,dwx,dwy,
                        dwx_b,dwy_b,pwx,pwy,pwx_b,pwy_b,theta,omega,alpha,temp_l};
int arr_len=sizeof(makeinit_array)/sizeof(double);

/*初期化ファイルの書き出し*/
initarray=fopen(name_initarray,"w");
for (i=0;i<N;i++){
    fprintf(initarray,"%lf    ",rx[i]);
    fprintf(initarray,"%lf    ",rx_b[i]);
    fprintf(initarray,"%lf    ",ry[i]);
    fprintf(initarray,"%lf    ",ry_b[i]);
    fprintf(initarray,"%lf    ",vx[i]);
    fprintf(initarray,"%lf    ",vy[i]);
    fprintf(initarray,"%lf    ",vy_b[i]);
    fprintf(initarray,"%lf    ",ax[i]);
    fprintf(initarray,"%lf    ",ay[i]);
    fprintf(initarray,"%lf    \n",ay_b[i]);
}
fclose(initarray);

initfile = fopen(name_init,"w");

for(i=0;i<arr_len;i++){
    if(i==arr_len-1) fprintf(initfile,"%lf",makeinit_array[i]);
    else fprintf(initfile,"%lf,",makeinit_array[i]);
}
fclose(initfile);
//file close
fclose(all_condition_file);
fclose(mini_file);
fclose(omega_file);
fclose(theta_file);
fclose(efile);
fclose(relation_heat_work);
return 0;
}

