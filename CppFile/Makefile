fwm_num.o:fwm_num.cpp
	g++ -o fwm_num.exe -c fwm_num.cpp

Funcs.o:Funcs.cpp
	g++ -o Funcs.exe -c Funcs.cpp

main: fwm_num.o Funcs.o
	g++ -o main fwm_num.exe Funcs.exe -lm -fopenmp