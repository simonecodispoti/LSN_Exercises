#include "funct.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
using namespace std;

int main(){
	unsigned int n_step = 10000;		//numero di step Montecarlo
	unsigned int n_cell = 100;		//numero di blocchi
	unsigned int l = n_step/n_cell;		//lunghezza dei blocchi

	Random ran;
	int* s = new int [4];
	
	for(int i=0; i<4; i++)
		s[i]=i;

	ran.SetRandom(s, 2, 7);			//parametri per LCG

	double* R = new double [n_step];

	for(int i=0; i<n_step; i++){		//array con numeri casuali distribuiti unif in [0,1)
		R[i] = ran.Rannyu();
	}

	double* ave = new double [n_cell];
	double* ave2 = new double [n_cell];
	double* sum_prog = new double [n_cell];		//progressione del valor medio
	double* sum2_prog = new double [n_cell];	//progressione della media dei quadrati
	double* err_prog = new double [n_cell];		//progressione dell'errore

	for(int i=0; i<n_cell; i++){
		double sum = 0;
			for(int j=0; j<l; j++){
				unsigned int pos = j+i*l;
				sum += R[pos];
			}
		ave[i] = sum/l;			//media per ogni cella
		ave2[i] = pow(ave[i],2);	//valore quadratico della media per ogni cella
	}

	for(int i=0; i<n_cell; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
	sum_prog[i]/=(i+1);
	sum2_prog[i]/=(i+1);

	if(i!=0)
		err_prog[i] = sqrt((sum2_prog[i] - pow(sum_prog[i],2))/i);	//incertezza statistica
	}

	err_prog[0] = 0;	//con un solo blocco non ho incertezza statistica

	Stampa("Mean.txt", sum_prog, n_cell);	//output: andamento del valor medio
	Stampa("Error.txt", err_prog, n_cell);	//output: andamento della deviazione standard della media

	return 0;
}














