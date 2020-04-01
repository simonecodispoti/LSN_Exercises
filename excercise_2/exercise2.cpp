#include "utilities.h"
#include "funzione.h"
#include "integraleMC.h"
using namespace std;

int main(){

	//  ||---|Stima di un integrale mediante metodi MC|---||  //

	int n_cell = 1000;		//Numero di calcoli dell'integrale
	int n_step = 1000;		//Numero di estrazioni per il calcolo di un singolo integrale	

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

	rnd.SaveSeed();

	//******************************RANDOM_GEN******************************//

	//---Metodo della media: Uniform sampling---//
    
	Coseno* cos = new Coseno(0.5*M_PI, 0.5*M_PI, 0);
	IntegraleMC integratore(rnd);

	double* esiti = new double [n_cell];

	for(int i=0; i<n_cell ;i++)
		esiti[i] = integratore.IntegraleMedia(0, 1, n_step, cos);
	
	Stampa("Integrali.txt", esiti, n_cell);

	double* media = new double [n_cell];
	MC_MeanProg(media, n_cell, n_cell, esiti);
	Stampa("I_unif.txt", media, n_cell);

	double* errore = new double [n_cell];
	MC_ErrProg(errore, n_cell, n_cell, esiti);
	Stampa("I_error_unif.txt", errore, n_cell);

	delete[] esiti;
	delete[] media;
	delete[] errore;


	//---Metodo della media: Importance sampling esponenziale---//

	Coseno_smorzato* f = new Coseno_smorzato(1-exp(-M_PI/2), 0.5*M_PI, 0.5*M_PI, 0);
	IntegraleMC integratore2;

	double* random = new double [n_step];
	double* esiti2 = new double [n_cell];

	for(int i=0; i<n_cell ;i++){
		for(int l=0; l<n_step; l++)
			random[l] = 0;
		for(int j=0; j<n_step; j++)
			random[j] = -(2/M_PI)*log(1-(1-exp(-M_PI/2))*rnd.Rannyu());
		esiti2[i] = integratore2.IntegraleMedia(0, 1, n_step, random, f);
	}
	
	Stampa("Integrali_exp.txt", esiti2, n_cell);

	double* media2 = new double [n_cell];
	MC_MeanProg(media2, n_cell, n_cell, esiti2);
	Stampa("I_exp.txt", media2, n_cell);

	double* errore2 = new double [n_cell];
	MC_ErrProg(errore2, n_cell, n_cell, esiti2);
	Stampa("I_error_exp.txt", errore2, n_cell);

	delete[] random;
	delete[] esiti2;
	delete[] media2;
	delete[] errore2;


	//  ||---|Random walks in 3D|---||  //


	return 0;

}
