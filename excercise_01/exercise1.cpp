#include "utilities.h"
#include "random.h"
using namespace std;

int main(){

	//  || --- | stima del valor medio e della varianza per una distribuzione uniforme | --- ||  //

	int n_step = 10000;		// numero di step Monte Carlo
	int n_cell = 100;		// numero di blocchi

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

	vector <double> R;

	for(int i=0; i<n_step; i++)		// array con numeri casuali distribuiti unif in [0,1)
		R.push_back(rnd.Rannyu());

	//******************************RANDOM_GEN******************************//

	// --- Media --- //

	vector <double> mean_prog;
	vector <double> err1_prog;

	MC_Mean_Error(R, mean_prog, err1_prog, n_step, n_cell);

	Print("Mean.txt", mean_prog);			// output: andamento del valor medio
	Print("Mean_error.txt", err1_prog);		// output: andamento dell'incertezza sulla media

	// --- Varianza --- //

	vector <double> V;

	for(int i=0; i<n_step; i++)				// la varaianza è la media del quadrato della distanza dal valor medio
		V.push_back(pow(R[i]-0.5,2));

	vector <double> var_prog;
	vector <double> err2_prog;

	MC_Mean_Error(V, var_prog, err2_prog, n_step, n_cell);

	Print("Var.txt", var_prog);				// output: andamento della varianza
	Print("Var_error.txt", err2_prog);		// output: andamento dell'incertezza sulla varianza


	//  || --- | test del Chi2 | --- ||  //

	// continuo ad usare n_step = 10000 e n_cell = 100; poichè il calcolo va effettuato 100 volte ho bisogno di n_step*n_cell numeri random

	int M = n_step*n_cell;
	int P_unif = n_step/n_cell;		// probabilità a priori che in n_step estrazioni un numero estratto si trovi in uno dei sottointervalli
	
	vector <double> R2;

	for(int i=0; i<M; i++)
		R2.push_back(rnd.Rannyu());	

	vector <double> atteso;
	vector <double> osservato(n_cell);
	vector <double> chi2(n_cell);

	for(int i=0; i<n_cell; i++)
		atteso.push_back(P_unif);
	
	for(int i=0; i<n_cell; i++){			// ciclo sul numero di simulazioni (n_cell)
		for(int l=0; l<n_cell; l++)			// ciclo di azzeramento 
			osservato[l] = 0;		
		for(int j=0; j<n_step; j++){		// ciclo interno ad ogni simulazione (n_step)
			int pos = j+i*n_step;
			for(int k=0; k<n_cell; k++){	// ciclo per valutare i conteggi in ogni "casella" che equipartiziona l'intervallo [0,1]
				if(R2[pos] > double(k)/n_cell && R2[pos] < double(k+1)/n_cell)
					osservato[k] += 1;
			}
		}
		chi2[i] = Chi2(osservato, atteso);
	}

	Print("Chi2.txt", chi2);


	//  || --- | Simulazione di un dado standard, esponenziale e lorentziano | --- ||  //

	// --- Prima simulazione: getto il dado N=10000 volte --- //

	int N = 10000;

	vector <double> G1;
	for(int i=0; i<N; i++)
		G1.push_back(rnd.Gauss(0,1));
	
	Print("Gauss_1.txt", G1);
	
	vector <double> E1;
	for(int i=0; i<N; i++)
		E1.push_back(rnd.Exp(1));

	Print("Exp_1.txt", E1);

	vector <double> L1;
	for(int i=0; i<N; i++)
		L1.push_back(rnd.Lorentz(0,1));

	Print("Lorentz_1.txt", L1);

	// --- Seconda simulazione: getto il dado N=2*10000=20000 volte e studio la variabile aleatoria x1 + x2 --- //

	vector <double> G2;
	for(int i=0; i<N; i++)
		G2.push_back((rnd.Gauss(0,1) + rnd.Gauss(0,1))/2);
	
	Print("Gauss_2.txt", G2);
	
	vector <double> E2;
	for(int i=0; i<N; i++)
		E2.push_back((rnd.Exp(1) + rnd.Exp(1))/2);

	Print("Exp_2.txt", E2);

	vector <double> L2;
	for(int i=0; i<N; i++)
		L2.push_back((rnd.Lorentz(0,1) + rnd.Lorentz(0,1))/2);

	Print("Lorentz_2.txt", L2);

	// --- Terza simulazione: getto il dado N=10*10000=10^5 volte e studio la variabile aleatoria x1 + x2 + ... + x10 --- //

	vector <double> G3;
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<10; j++)
			sum += rnd.Gauss(0,1);
		sum/=10;
		G3.push_back(sum);
	}

	Print("Gauss_3.txt", G3);
	
	vector <double> E3;
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<10; j++)
			sum += rnd.Exp(1);
		sum/=10;
		E3.push_back(sum);
	}

	Print("Exp_3.txt", E3);

	vector <double> L3;
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<10; j++)
			sum += rnd.Lorentz(0,1);
		sum/=10;
		L3.push_back(sum);
	}

	Print("Lorentz_3.txt", L3);

	// --- Quarta simulazione: getto il dado N=100*10000=10^6 volte e studio la variabile aleatoria x1 + x2 + ... +x100 --- //

	vector <double> G4;
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<100; j++)
			sum += rnd.Gauss(0,1);
		sum/=100;
		G4.push_back(sum);
	}

	Print("Gauss_4.txt", G4);

	vector <double> E4;
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<100; j++)
			sum += rnd.Exp(1);
		sum/=100;
		E4.push_back(sum);
	}

	Print("Exp_4.txt", E4);

	vector <double> L4;
	for(int i=0; i<N; i++){
		double sum = 0;
		for(int j=0; j<100; j++)
			sum += rnd.Lorentz(0,1);
		sum/=100;
		L4.push_back(sum);
	}

	Print("Lorentz_4.txt", L4);


	//  || --- | Simulazione dell'esperimento di Buffon | --- ||  //

	// supponiamo di lanciare l'ago n_step = 1000 volte e ripetiamo l'esperimento n_cell = 100 volte ---> effettuo in totale n_step*n_cell = 100000 lanci 

	n_step = 1000;
	n_cell = 100;

	double l = 4.5;		//lunghezza dell'ago in cm
	double d = 7.0;		//distanza tra gli assi del parquet in cm :)

	vector <double> Pi;
	
	int n_hit = 0;

	// sqrt(pow(l,2)/4 - pow(rnd.Rannyu(0,l/2),2))) non funziona!?

	for(int i=0; i<n_cell; i++){
		n_hit = 0;
		for(int j=0; j<n_step; j++)
			if (rnd.Rannyu(0,d/2) <= l/2*sin(rnd.Rannyu(0,0.5*M_PI))) n_hit ++;
	Pi.push_back((2*l*n_step)/(d*n_hit));
	}

	vector <double> media;
	vector <double> errore;

	MC_Mean_Error(Pi, media, errore, n_cell, n_cell);

	Print("Pi.txt", media);
	Print("Pi_error.txt", errore);

	rnd.SaveSeed();		// salvataggio del seme per riproducibilità

	return 0;
}
