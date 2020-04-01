#include "utilities.h"
#include "random.h"
using namespace std;

int main(){

	//  ||---|stima del valor medio e della varianza per una distribuzione uniforme|---||  //

	int n_step = 10000;		//numero di step Monte Carlo
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


	//  ||---|test del Chi2|---||  //

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
	
	for(int i=0; i<n_cell; i++){			//ciclo sul numero di simulazioni (n_cell)
		for(int l=0; l<n_cell; l++)		
			osservato[l] = 0;		
		for(int j=0; j<n_step; j++){		//ciclo interno ad ogni simulazione (n_step)
			int pos = j+i*n_step;
			for(int k=0; k<n_cell; k++){	//ciclo per valutare i conteggi in ogni "casella" che equipartiziona l'intervallo [0,1]
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


	//  ||---|Simulazione di un dado standard, esponenziale e lorentziano|---||  //

	//Prima simulazione: getto il dado N=10000 volte

	int N = 10000;

	double* G1 = new double [N];
	for(int i=0; i<N; i++)
		G1[i] = rnd.Gauss(0,1);
	
	Stampa("Gauss_1.txt", G1, N);
	delete[] G1;
	
	double* E1 = new double [N];
	for(int i=0; i<N; i++)
		E1[i] = rnd.Exp(1);

	Stampa("Exp_1.txt", E1, N);
	delete[] E1;

	double* L1 = new double [N];
	for(int i=0; i<N; i++)
		L1[i] = rnd.Lorentz(0,1);

	Stampa("Lorentz_1.txt", L1, N);
	delete[] L1;

	//Seconda simulazione: getto il dado N=2*10000=20000 volte e studio la variabile aleatoria x1 + x2

	double* G2 = new double [N];
	for(int i=0; i<N; i++)
		G2[i] = (rnd.Gauss(0,1) + rnd.Gauss(0,1))/2;
	
	Stampa("Gauss_2.txt", G2, N);
	delete[] G2;
	
	double* E2 = new double [N];
	for(int i=0; i<N; i++)
		E2[i] = (rnd.Exp(1) + rnd.Exp(1))/2;

	Stampa("Exp_2.txt", E2, N);
	delete[] E2;

	double* L2 = new double [N];
	for(int i=0; i<N; i++)
		L2[i] = (rnd.Lorentz(0,1) + rnd.Lorentz(0,1))/2;

	Stampa("Lorentz_2.txt", L2, N);
	delete[] L2;

	//Terza simulazione: getto il dado N=10*10000=10^5 volte e studio la variabile aleatoria x1 + x2 + ... + x10

	double* G3 = new double [N];
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<10; j++)
			sum += rnd.Gauss(0,1);
		sum/=10;
		G3[i] = sum;
	}

	Stampa("Gauss_3.txt", G3, N);
	delete[] G3;
	
	double* E3 = new double [N];
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<10; j++)
			sum += rnd.Exp(1);
		sum/=10;
		E3[i] = sum;
	}

	Stampa("Exp_3.txt", E3, N);
	delete[] E3;

	double* L3 = new double [N];
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<10; j++)
			sum += rnd.Lorentz(0,1);
		sum/=10;
		L3[i] = sum;
	}

	Stampa("Lorentz_3.txt", L3, N);
	delete[] L3;

	//Quarta simulazione: getto il dado N=100*10000=10^6 volte e studio la variabile aleatoria x1 + x2 + ... +x100

	double* G4 = new double [N];
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<100; j++)
			sum += rnd.Gauss(0,1);
		sum/=100;
		G4[i] = sum;
	}

	Stampa("Gauss_4.txt", G4, N);
	delete[] G4;

	double* E4 = new double [N];
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<100; j++)
			sum += rnd.Exp(1);
		sum/=100;
		E4[i] = sum;
	}

	Stampa("Exp_4.txt", E4, N);
	delete[] E4;

	double* L4 = new double [N];
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<100; j++)
			sum += rnd.Lorentz(0,1);
		sum/=100;
		L4[i] = sum;
	}

	Stampa("Lorentz_4.txt", L4, N);
	delete[] L4;


	//  ||---|Simulazione dell'esperimento di Buffon|---||  //

	//supponiamo di lanciare l'ago n_step = 1000 volte e ripetiamo l'esperimento n_cell = 100 volte ---> effettuo in totale n_step*n_cell = 100000 lanci 

	n_step = 1000;
	n_cell = 100;

	double l = 4.5;		//lunghezza dell'ago in cm
	double d = 7.0;		//distanza tra gli assi del parquet in cm :)

	double* Pi = new double [n_cell];
	
	int n_hit = 0;

	//sqrt(pow(l,2)/4 - pow(rnd.Rannyu(0,l/2),2))) non funziona?

	for(int i=0; i<n_cell; i++){
		n_hit = 0;
		for(int j=0; j<n_step; j++){
			if (rnd.Rannyu(0,d/2) <= l/2*sin(rnd.Rannyu(0,0.5*M_PI)))
				n_hit++;
		}
	Pi[i] = (2*l*n_step)/(d*n_hit);
	}

	double* media = new double [n_cell];
	MC_MeanProg(media, n_cell, n_cell, Pi);
	Stampa("Pi.txt", media, n_cell);

	double* errore = new double [n_cell];
	MC_ErrProg(errore, n_cell, n_cell, Pi);
	Stampa("Pi_error.txt", errore, n_cell);

	delete[] Pi;
	delete[] media;
	delete[] errore;

	return 0;
}
