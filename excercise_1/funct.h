#ifndef FUNCT_H
#define FUNCT_H
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
using namespace std;

void Stampa(const char* filename, const double* v, const int dim);				//stampa double su file
void MC_MeanProg(double* media, const int n_step, const int n_cell, const double* data);	//Andamento della miglior stima Montecarlo di una grandezza double in una simulazione di n_cell lunghe n_step
void MC_ErrProg(double* errore, const int n_step, const int n_cell, const double* data);	//Andamento dell'errore per la miglior stima Montecarlo ... 
double Chi2(const double* obs, const double* exp, const int dim);				//test del Chi2 (dim è la dimensione del vettore delle frequenze attese, che deve concidere con quella delle 													frequenze osservate)

#endif /* FUNCT_H */
