#include "funct.h"
#include "random.h"
using namespace std;

int main(){

	//---|stima del valor medio e della varianza per una distribuzione uniforme|---//

	int n_step = 10000;		//numero di step Montecarlo
	int n_cell = 100;		//numero di blocchi

	//******************************RANDOM_GEN******************************//
	
	Random rnd;
	int seed[4];
	int p1, p2;

	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
				if( property == "RANDOMSEED" ){
					input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
					rnd.SetRandom(seed,p1,p2);
				}
		}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	double* R = new double [n_step];

	for(int i=0; i<n_step; i++)		//array con numeri casuali distribuiti unif in [0,1)
		R[i] = rnd.Rannyu();

	rnd.SaveSeed();

	//******************************RANDOM_GEN******************************//

	//---Media---//

	double* mean_prog = new double [n_cell];
	MC_MeanProg(mean_prog, n_step, n_cell, R);
	Stampa("Mean.txt", mean_prog, n_cell);			//output: andamento del valor medio
	
	double* err1_prog = new double [n_cell];
	MC_ErrProg(err1_prog, n_step, n_cell, R);
	Stampa("Mean_error.txt", err1_prog, n_cell);		//output: andamento dell'incertezza sulla media

	delete[] mean_prog;
	delete[] err1_prog;

	//---Varianza---//

	double* V = new double [n_step];

	for(int i=0; i<n_step; i++)				//la varaianza è la media del quadrato della distanza dall valor medio
		V[i] = pow(R[i]-0.5,2);

	double* var_prog = new double [n_cell];
	MC_MeanProg(var_prog, n_step, n_cell, V);
	Stampa("Var.txt", var_prog, n_cell);			//output: andamento della varianza

	double* err2_prog = new double [n_cell];
	MC_ErrProg(err2_prog, n_step, n_cell, V);
	Stampa("Var_error.txt", err2_prog, n_cell);		//output: andamento dell'incertezza sulla varianza

	delete[] R;
	delete[] V;
	delete[] var_prog;
	delete[] err2_prog;


	//---|test del Chi2|---//

	//continuo ad usare n_step = 10000 e n_cell = 100; poichè il calcolo va effettuato 100 volte ho bisogno di n_step*n_cell numeri random

	int M = n_step*n_cell;
	int P_unif = n_step/n_cell;	//probabilità a priori che in n_step estrazioni un numero estratto si trovi in uno dei sottointervalli
	
	double* R2 = new double [M];

	for(int i=0; i<M; i++)
		R2[i] = rnd.Rannyu();	

	double* atteso = new double [n_cell];
	double* osservato = new double [n_cell];
	double* chi2 = new double [n_cell];

	for(int i=0; i<n_cell; i++)
		atteso[i] = P_unif;
	
	for(int i=0; i<n_cell; i++){
		for(int l=0; l<n_cell; l++)
			osservato[l] = 0;
		for(int j=0; j<n_step; j++){
			int pos = j+i*n_step;
			for(int k=0; k<n_cell; k++){
				if(R2[pos] > double(k)/n_cell && R2[pos] < double(k+1)/n_cell)
					osservato[k] += 1;
			}
		}
		chi2[i] = Chi2(osservato, atteso, n_cell);
	}

	Stampa("Chi2.txt", chi2, n_cell);
	
	delete[] R2;
	delete[] atteso;
	delete[] osservato; 
	delete[] chi2;

	return 0;
}












