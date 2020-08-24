CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise9.exe : exercise9.o random.o posizione.o TSP.o utilities.h
	$(CC) exercise9.o random.o posizione.o TSP.o -o exercise9.exe

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

posizione.o: posizione.cpp posizione.h
		$(CC) -c posizione.cpp -o posizione.o $(CFLAGS)

TSP.o : TSP.cpp TSP.h 
	$(CC) -c TSP.cpp -o TSP.o $(CFLAGS)

clean :
	rm *.o exercise9.exe