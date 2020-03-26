CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise1.exe : exercise1.o random.o funct.o
	$(CC) exercise1.o random.o funct.o -o exercise1.exe

funct.o : funct.cpp funct.h
	$(CC) -c funct.cpp -o funct.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o exercise1.exe Mean.txt Mean_error.txt Var.txt Var_error.txt Chi2.txt
