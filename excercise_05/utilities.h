#ifndef UTILITIES_H
#define UTILITIES_H
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
using namespace std;

void Stampa(const char* filename, const double* v, const int dim);				//Stampa vettore di double su file

void MC_MeanProg(double* media, const int n_step, const int n_cell, const double* data);	//Andamento della miglior stima Monte Carlo di una grandezza in una simulazione di n_cell blocchi di 													lunghezza pari a n_step; bisogna sempre interpretare n_step MC in base alla specificità della simulazione: 													quasi mai coniciderà con la singola estrazione di un numero casuale, e ciò porterà a porre n_step = n_cell

void MC_ErrProg(double* errore, const int n_step, const int n_cell, const double* data);	//Andamento dell'errore relativo al metodo MC_MeanProg

void Eval_ave_err(const double* input, double* average, double* error, const int n_step, const int n_cell);	//metodo misto: posso valutare media ed errore simultaneamente

double Chi2(const double* obs, const double* exp, const int dim);				//Test del Chi2 (dim è la dimensione del vettore delle frequenze attese, che deve concidere con quella delle 													frequenze osservate)

#endif /* UTILITIES_H */
