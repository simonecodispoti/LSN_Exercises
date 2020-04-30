CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise5.exe : exercise5.o random.o utilities.o posizione.o funzione.o integraleMC.o
	$(CC) exercise5.o random.o utilities.o posizione.o funzione.o integraleMC.o -o exercise5.exe

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

utilities.o : utilities.cpp utilities.h
	$(CC) -c utilities.cpp -o utilities.o $(CFLAGS)

posizione.o : posizione.cpp posizione.h
	$(CC) -c posizione.cpp -o posizione.o $(CFLAGS)

funzione.o : funzione.cpp funzione.h
	$(CC) -c funzione.cpp -o funzione.o $(CFLAGS)

integraleMC.o : integraleMC.cpp integraleMC.h
	$(CC) -c integraleMC.cpp -o integraleMC.o $(CFLAGS)

clean :
	rm *.o exercise5.exe
