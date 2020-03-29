#ifndef UTILITIES_H
#define UTILITIES_H
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
using namespace std;

void Stampa(const char* filename, const double* v, const int dim);				//stampa vettore di double su file

void MC_MeanProg(double* media, const int n_step, const int n_cell, const double* data);	//Andamento della miglior stima Monte Carlo di una grandezza in una simulazione di n_cell blocchi di 													lunghezza pari a n_step; se si deve stimare l'andamento di una grandezza per la quale si sono già effettuati 													n_step attraverso metodi esterni (p.e. valutazione di un integrale col metodo della media) basta porre 													n_step == n_cell

void MC_ErrProg(double* errore, const int n_step, const int n_cell, const double* data);	//Andamento dell'errore relativo al metodo MC_MeanProg

double Chi2(const double* obs, const double* exp, const int dim);				//test del Chi2 (dim è la dimensione del vettore delle frequenze attese, che deve concidere con quella delle 													frequenze osservate)

#endif /* UTILITIES_H */
