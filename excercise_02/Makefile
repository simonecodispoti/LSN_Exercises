CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise2.exe : exercise2.o random.o funzione.o integraleMC.o posizione.o utilities.h
	$(CC) exercise2.o random.o funzione.o integraleMC.o posizione.o -o exercise2.exe

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

funzione.o : funzione.cpp funzione.h
	$(CC) -c funzione.cpp -o funzione.o $(CFLAGS)

integraleMC.o : integraleMC.cpp integraleMC.h
	$(CC) -c integraleMC.cpp -o integraleMC.o $(CFLAGS)

posizione.o: posizione.cpp posizione.h
	$(CC) -c posizione.cpp -o posizione.o $(CFLAGS)

clean :
	rm *.o exercise2.exe Integrali.txt I_unif.txt I_error_unif.txt Integrali_exp.txt I_exp.txt I_error_exp.txt dist_lattice.txt dist_err_lattice.txt dist_continuum.txt dist_err_continuum.txt X.txt Y.txt Z.txt Xc.txt Yc.txt Zc.txt seed.out
