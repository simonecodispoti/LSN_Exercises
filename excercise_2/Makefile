CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise2.exe : exercise2.o random.o utilities.o funzione.o integraleMC.o
	$(CC) exercise2.o random.o utilities.o funzione.o integraleMC.o -o exercise2.exe

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

utilities.o : utilities.cpp utilities.h
	$(CC) -c utilities.cpp -o utilities.o $(CFLAGS)

funzione.o : funzione.cpp funzione.h
	$(CC) -c funzione.cpp -o funzione.o $(CFLAGS)

integraleMC.o : integraleMC.cpp integraleMC.h
	$(CC) -c integraleMC.cpp -o integraleMC.o $(CFLAGS)

clean :
	rm *.o exercise2.exe Integrali.txt I_unif.txt I_error_unif.txt Integrali_exp.txt I_exp.txt I_error_exp.txt
