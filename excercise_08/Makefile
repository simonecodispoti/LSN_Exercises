CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise8.exe : exercise8.o random.o VQMC.o utilities.h
	$(CC) exercise8.o random.o VQMC.o -o exercise8.exe

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

VQMC.o : VQMC.cpp VQMC.h
	$(CC) -c VQMC.cpp -o VQMC.o $(CFLAGS)

clean :
	rm *.o exercise8.exe seed.out