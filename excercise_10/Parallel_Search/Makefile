CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

exercise10_Parallel.exe : exercise10_Parallel.o random.o posizione.o TSP.o utilities.h
	$(CC) exercise10_Parallel.o random.o posizione.o TSP.o -o exercise10_Parallel.exe

exercise10_Parallel.o : exercise10_Parallel.cpp
	$(CC) -c exercise10_Parallel.cpp -o exercise10_Parallel.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

posizione.o: posizione.cpp posizione.h
	$(CC) -c posizione.cpp -o posizione.o $(CFLAGS)

TSP.o : TSP.cpp TSP.h 
	$(CC) -c TSP.cpp -o TSP.o $(CFLAGS)

clean :
	rm *.o exercise10_Parallel.exe seed.out
