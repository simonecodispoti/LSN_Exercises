CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise3.exe : exercise3.o random.o utilities.o
	$(CC) exercise3.o random.o utilities.o -o exercise3.exe

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

utilities.o : utilities.cpp utilities.h
	$(CC) -c utilities.cpp -o utilities.o $(CFLAGS)

clean :
	rm *.o exercise3.exe Call_direct.txt Call_direct_err.txt Put_direct.txt Put_direct_err.txt Call_discret.txt Call_discret_err.txt Put_discret.txt Put_discret_err.txt
