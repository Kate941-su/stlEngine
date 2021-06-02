#include "modelC_init.h"
#define MT_N 624
#define MT_M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[MT_N]; /* the array for the state vector  */
static int mti = MT_N + 1; /* mti==MT_N+1 means mt[MT_N] is not initialized */

///////////メルセンヌツイスター///////////
/* initializes mt[MT_N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0] = s & 0xffffffffUL;
    for (mti = 1; mti < MT_N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i = 1; j = 0;
    k = (MT_N > key_length ? MT_N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL))
            + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
        if (j >= key_length) j = 0;
    }
    for (k = MT_N - 1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL))
            - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i >= MT_N) { mt[0] = mt[MT_N - 1]; i = 1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2] = { 0x0UL, MATRIX_A };
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= MT_N) { /* generate N words at one time */
        int kk;

        if (mti == MT_N + 1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk = 0;kk < MT_N - MT_M;kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk < MT_N - 1;kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MT_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32() * (1.0 / 4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32() * (1.0 / 4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5) * (1.0 / 4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
    unsigned long a = genrand_int32() >> 5, b = genrand_int32() >> 6;
    return(a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */
//////////////////////////////////////////////////////////////////


//doubleの最大値を返す
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

//int型の合計値を返す
int int_array_sum(int array[N]){
    int i;
    int sum=0;
    for (i=0;i<N;i++){
        sum+=array[i];
    }
    return sum;
}

//double型の合計値を返す
double double_array_sum(double array[N]){
    int i;
    double sum=0;
    for (i=0;i<N;i++){
        sum+=array[i];
    }
    return sum;
}

//不明
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

//周期境界条件（左右境界）
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

//弾性壁（上下）
void elastic(int NP, double RY[], double VY[], double LY) {
	int i;
	for (i = 0;i < NP;i++) {
		if (RY[i] > LY - 0.5 | RY[i] < 0.5) {
			VY[i] = -VY[i];
		}
	}
}

//グリッドを系に設置する
void gmap_create(int NP, double RX[], double RY[], double L_GX, double L_GY,
                 int N_GX, int N_GY, int NEIGHBOR_ROW[], int NEIGHBOR_COL[], 
                 int NEITGHBOR_LEN, int PAIRLIST[][10], double LX, int **G_MAP) {
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

//力の計算をする(高速化無し)
void force_nos(int NP, double RX[], double RY[], double AX[], double AY[], double LX, 
                double LY, double RC2, double* POT, double* POT_IJ, double* MINI) {
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

//力の計算をする
void force(int NP, double RX[], double RY[], double AX[], double AY[], double LX, 
            double LY, int PAIRLIST[][10], double RC2, double* POT, 
            double* POT_IJ, double* MINI) {
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

//熱壁（下固定、上パワーピストン）
void heatwall(double H,int NP,double RY[],double RY_B[],double VY[],double *Q_IN,
              double *Q_OUT ,double PPY,double PPV,double TEMP_L,double TEMP_H,
              double *FPP,double LY,double *h_ss,double *d_w,double *q_debug){
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

//上から下へ粒子が通過するときに熱の計算を行う
void piston_move_u(int NP,double RY[N],double RY_B[N],double VY[N],double VY_B[N],
                    double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,
                    double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,
                    int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,
                    double *H1_D,int *k,int *j,double MDP,double PROBABIRITY,double *q_debug){
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

//下から上へ粒子が通過するときに熱の計算を行う
void piston_move_d(int NP,double RY[N],double RY_B[N],double VY[N],double VY_B[N],
                    double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,
                    double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,
                    int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,
                    double *H1_D,int *kk,int *jj,double MDP,double PROBABIRITY,double *q_debug){
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

//温度の計算を行う
void calculateTemleture(int NP,double RY[],double VX[],double VY[],
                        double DPY,double *TEMP_D,double *TEMP_U,int *ND,int*NU ,
                        int N_U_list[N],double *TEMP){
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

//温度と圧力を測る
void calcTempAndPress(double RY[],double VX[] , double VY[] ,double PPY,
                double TEMP_PART_LIST[NWRITE][PARTNUM],
                double PRESS_PART_LIST[NWRITE][PARTNUM],int T_WRITE){
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