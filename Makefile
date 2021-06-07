#初期化ファイル作成
#init_N : ヘッダファイルで粒子数を操作する
init_N:
	gcc -o modelCinitN.exe modelC_init_N.c modelC_init_sub.c -lm -fopenmp
#init_temp : コマンドライン引数で温度を操作する
init_temp:
	gcc -o modelCinittemp.exe modelC_init_temp.c modelC_init_sub.c -lm -fopenmp
#init_prob : コマンドライン引数で確率を操作する
init_prob:
	gcc -o modelCinittemp.exe modelC_init_temp.c modelC_init_sub.c -lm -fopenmp

#定常ファイル作成
#stable_N : ヘッダファイルで粒子数を操作する
stable_N:
	gcc -o modelCstableN.exe modelC_stable_N.c modelC_init_sub.c -lm -fopenmp
#init_temp : コマンドライン引数で温度を操作する
stable_temp:
	gcc -o modelCstabeltemp.exe modelC_stable_temp.c modelC_init_sub.c -lm -fopenmp
#init_temp : コマンドライン引数で温度を操作する
stable_prob:
	gcc -o modelCstableprob.exe modelC_stable_prob.c modelC_init_sub.c -lm -fopenmp
