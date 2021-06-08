#pragma once
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define Nx 30
#define Ny 15
#define PARTNUM 2
#define N (Nx*Ny)
//server
#define NSTEPS 1000000
//local test
//#define NSTEPS 200000
#define NDEVIDE 100000
#define NWRITE (NSTEPS/NDEVIDE)

//////メルセンヌツイスタ///////
void init_genrand(unsigned long s);

void init_by_array(unsigned long init_key[], int key_length);

unsigned long genrand_int32(void);

long genrand_int31(void);

double genrand_real1(void);

double genrand_real2(void);

double genrand_real3(void);

double genrand_res53(void);
/////////////////////////////////

//doubleの最大値を返す
double judge_max(double array[],int len_array);

//int型の合計値を返す
int int_array_sum(int array[N]);

//double型の合計値を返す
double double_array_sum(double array[N]);

//不明
int bigger(int a,int b,int *UP_DOWN_FLAG);

//周期境界条件（左右境界）
void boundary(int NP, double RX[], double LX);

//弾性壁（上下）
void elastic(int NP, double RY[], double VY[], double LY);

//グリッドを系に設置する
void gmap_create(int NP, double RX[], double RY[], double L_GX, double L_GY, int N_GX, 
                int N_GY, int NEIGHBOR_ROW[], int NEIGHBOR_COL[], int NEITGHBOR_LEN, 
                int PAIRLIST[][10], double LX, int **G_MAP);

//力の計算をする
void force(int NP, double RX[], double RY[], double AX[], double AY[], double LX, 
            double LY, int PAIRLIST[][10], double RC2, double* POT, double* POT_IJ, double* MINI);

//力の計算をする(高速化無し)
void force_nos(int NP, double RX[], double RY[], double AX[], double AY[], double LX,
                 double LY, double RC2, double* POT, double* POT_IJ, double* MINI);

//熱壁（下固定、上パワーピストン）
void heatwall(double H,int NP,double RY[],double RY_B[],double VY[],double *Q_IN,
                double *Q_OUT ,double PPY,double PPV,double TEMP_L,double TEMP_H,
                double *FPP,double LY,double *h_ss,double *d_w);

//上から下へ粒子が通過するときに熱の計算を行う
void piston_move_u(int NP,double RY[N],double RY_B[N],double VY[N],double VY_B[N],
                    double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,
                    double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,
                    int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,
                    double *H1_D,int *k,int *j,double MDP,double PROBABIRITY);

//下から上へ粒子が通過するときに熱の計算を行う
void piston_move_d(int NP,double RY[N],double RY_B[N],double VY[N],double VY_B[N],
                    double AY[],double H,double H_REV,double *Q_IN,double *Q_OUT,
                    double Q_IN_SUM,double Q_OUT_SUM,double DPY,double DPY_B,double *DPV,
                    int *HIT_PISTON,int *THROUGH_PISTON,double *FDP,double TEMP_L,double TEMP_H,
                    double *H1_D,int *kk,int *jj,double MDP,double PROBABIRITY);

//温度の計算を行う
void calculateTemleture(int NP,double RY[],double VX[],double VY[],double DPY,double *TEMP_D,
                double *TEMP_U,int *ND,int*NU ,int N_U_list[N],double *TEMP);

//温度と圧力を測る
void calcTempAndPress(double RY[],double VX[] , double VY[] ,double PPY,double TEMP_PART_LIST[NWRITE][PARTNUM],
                double PRESS_PART_LIST[NWRITE][PARTNUM],int T_WRITE);
