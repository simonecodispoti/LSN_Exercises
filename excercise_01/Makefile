CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise1.exe : exercise1.o random.o utilities.h
	$(CC) exercise1.o random.o -o exercise1.exe

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o exercise1.exe seed.out Mean.txt Mean_error.txt Var.txt Var_error.txt Chi2.txt Gauss_1.txt Gauss_2.txt Gauss_3.txt Gauss_4.txt Exp_1.txt Exp_2.txt Exp_3.txt Exp_4.txt Lorentz_1.txt Lorentz_2.txt Lorentz_3.txt Lorentz_4.txt Pi.txt Pi_error.txt seed.out