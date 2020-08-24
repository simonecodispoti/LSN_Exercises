CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise5.exe : exercise5.o random.o posizione.o funzione.o utilities.h
	$(CC) exercise5.o random.o posizione.o funzione.o -o exercise5.exe

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

posizione.o : posizione.cpp posizione.h
	$(CC) -c posizione.cpp -o posizione.o $(CFLAGS)

funzione.o : funzione.cpp funzione.h
	$(CC) -c funzione.cpp -o funzione.o $(CFLAGS)

clean :
	rm *.o exercise5.exe X_s.txt X_p.txt Y_s.txt Y_p.txt Z_s.txt Z_p.txt Mean_unif_s.txt Mean_unif_p.txt Mean_norm_s.txt Mean_norm_p.txt Error_unif_s.txt Error_unif_p.txt Error_norm_s.txt Error_norm_p.txt seed.out
