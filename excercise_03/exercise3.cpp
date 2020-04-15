#include "random.h"
#include "utilities.h"
using namespace std;

int main(){

	//  ||---|Stima del prezzo di un'opzione europea|---||  //

	int n_cell = 100;	//Numero di simulazioni
	int n_step = 10000;	//Numero di step Monte Carlo: è la simulazione dell'intera traiettoria di un GBM

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

	//******************************RANDOM_GEN******************************//
	
	//Parametri del contratto

	double S0 = 100;
	double K = 100;
	double T = 1;

	//Parametri di mercato

	double r = 0.1;
	double sigma = 0.25;


	//---Calcolo diretto delle opzioni Call e Put: n_step GBM(r,sigma^2) tramite salto diretto da 0 a T---//

	double* call = new double [n_step];
	double* put = new double [n_step];
	double S = 0;
	
	for(int i=0; i<n_step; i++){
		S = S0*exp((r-0.5*pow(sigma,2))*T + sigma*rnd.Gauss(0,1)*sqrt(T));
		if(S > K){
			call [i] = exp(-r*T)*(S-K);
			put [i] = 0;
		}
		if(S < K){
			call [i] = 0;
			put [i] = exp(-r*T)*(K-S);	
		}
	}
	
	double* call_media = new double [n_cell];
	MC_MeanProg(call_media, n_step, n_cell, call);
	Stampa("Call_direct.txt", call_media, n_cell);

	double* call_error = new double [n_cell];
	MC_ErrProg(call_error, n_step, n_cell, call);
	Stampa("Call_direct_err.txt", call_error, n_cell);

	double* put_media = new double [n_cell];
	MC_MeanProg(put_media, n_step, n_cell, put);
	Stampa("Put_direct.txt", put_media, n_cell);

	double* put_error = new double [n_cell];
	MC_ErrProg(put_error, n_step, n_cell, put);
	Stampa("Put_direct_err.txt", put_error, n_cell);


	//---Calcolo con discretizzazione delle opzioni Call e Put: n_step GBM(0,sigma^2) tramite passi successivi di ampiezza dt da 0 a T---//
	
	for(int l=0; l<n_step; l++){
		call [l] = 0;
		put [l] = 0;
	}

	for(int l=0; l<n_cell; l++){
		call_media [l] = 0;
		call_error [l] = 0;
		put_media [l] = 0;
		put_error [l] = 0;
	}

	int L = 100;		//numero di intervalli temporali in [0,T]
	double dt = T/L;	//passo temporale
	double S_prec = S0;	//variabile d'appoggio per memorizzare il valore di S(t)
	S = 0;
	
	for(int i=0; i<n_step; i++){
		S_prec = S0;
		for(int j=0; j<L; j++){
			S = S_prec*exp((r-0.5*pow(sigma,2))*dt + sigma*rnd.Gauss(0,1)*sqrt(dt));
			S_prec = S;
		}
		if(S > K){
			call [i] = exp(-r*T)*(S-K);
			put [i] = 0;
		}
		if(S < K){
			call [i] = 0;
			put [i] = exp(-r*T)*(K-S);	
		}
	}

	MC_MeanProg(call_media, n_step, n_cell, call);
	Stampa("Call_discret.txt", call_media, n_cell);

	MC_ErrProg(call_error, n_step, n_cell, call);
	Stampa("Call_discret_err.txt", call_error, n_cell);

	MC_MeanProg(put_media, n_step, n_cell, put);
	Stampa("Put_discret.txt", put_media, n_cell);

	MC_ErrProg(put_error, n_step, n_cell, put);
	Stampa("Put_discret_err.txt", put_error, n_cell);

	delete[] call;
	delete[] put;
	delete[] call_media;
	delete[] call_error;
	delete[] put_media;
	delete[] put_error;

	rnd.SaveSeed();		//Salvo il seme per riproducibilità

	return 0;
}
