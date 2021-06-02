hello: modelC_init_temp.c modelC_init.c modelC_init.h
	gcc -o sample.exe modelC_init_temp.c modelC_init.c -lm -fopenmp