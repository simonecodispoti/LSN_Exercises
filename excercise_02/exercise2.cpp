#include "utilities.h"
#include "funzione.h"
#include "integraleMC.h"
#include "posizione.h"
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

	//******************************RANDOM_GEN******************************//

	//---Metodo della media: Uniform sampling---//
    
	Coseno* Cos = new Coseno(0.5*M_PI, 0.5*M_PI, 0);
	IntegraleMC integratore(rnd);

	double* esiti = new double [n_cell];

	for(int i=0; i<n_cell ;i++)
		esiti[i] = integratore.IntegraleMedia(0, 1, n_step, Cos);
	
	Stampa("Integrali.txt", esiti, n_cell);

	double* media = new double [n_cell];
	MC_MeanProg(media, n_cell, n_cell, esiti);
	Stampa("I_unif.txt", media, n_cell);

	double* errore = new double [n_cell];
	MC_ErrProg(errore, n_cell, n_cell, esiti);
	Stampa("I_error_unif.txt", errore, n_cell);

	delete Cos;
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

	delete f;
	delete[] random;
	delete[] esiti2;
	delete[] media2;
	delete[] errore2;


	//  ||---|Random walks in 3D|---||  //

	//---Reticolo di spaziatura unitaria---//

	double a = 1;		//spaziatura del reticolo
	n_step = 1000;		//numero di passi del camminatore
	n_cell = 10000;		//Numero di RW

	Posizione O;		//origine (0,0,0)

	Posizione** RW = new Posizione* [n_cell];		//vettore di RW lungo n_cell
	
	for(int i=0; i<n_cell; i++)				//vettori di Posizioni lunghi n_step
		RW[i] = new Posizione [n_step];

	for(int i=0; i<n_cell; i++){				//inizializzio tutte le posizioni a (0,0,0)
		for(int j=0; j<n_step; j++)
			RW[i][j] = O;
	}

	double p = 0;						//variabili di appoggio
	double sum_x = 0;
	double sum_y = 0;
	double sum_z = 0;

	for(int i=0; i<n_cell; i++){				//ciclo sulle n_cell RW
		sum_x = 0;
		sum_y = 0;
		sum_z = 0;
		for(int j=1; j<n_step; j++){			//ciclo sugli n_step di ogni RW
			p = rnd.Rannyu();
			if(p < 0.5){				//campionamento di variabile aleatoria discreta tramite partizione dell'intervallo [0,1)
				if(p < double(1)/3){
					if(p < double(1)/6){
						sum_x += a;
						RW[i][j].SetX(sum_x);
						RW[i][j].SetY(sum_y);
						RW[i][j].SetZ(sum_z);
						continue;
					}
				sum_x -= a;
				RW[i][j].SetX(sum_x);
				RW[i][j].SetY(sum_y);
				RW[i][j].SetZ(sum_z);
				continue;
				}
			sum_y += a;
			RW[i][j].SetX(sum_x);
			RW[i][j].SetY(sum_y);
			RW[i][j].SetZ(sum_z);
			}
			else if (p > 0.5){
				if(p > double(2)/3){
					if(p > double(5)/6){
						sum_z -= a;
						RW[i][j].SetX(sum_x);
						RW[i][j].SetY(sum_y);
						RW[i][j].SetZ(sum_z);
						continue;
					}
				sum_z += a;
				RW[i][j].SetX(sum_x);
				RW[i][j].SetY(sum_y);
				RW[i][j].SetZ(sum_z);
				continue;
				}
			sum_y -= a;
			RW[i][j].SetX(sum_x);
			RW[i][j].SetY(sum_y);
			RW[i][j].SetZ(sum_z);
			}
		}
	}

	int N = n_step*n_cell;
	double* dist_quad = new double [N];		//vettore delle distanze quadratiche dall'origine

	for(int i=0; i<n_step; i++){					//Distanza allo step i-esimo...
		for(int j=0; j<n_cell; j++)				//...nella j-esima RW
			dist_quad[i*n_cell + j] = RW[j][i].Distanza_quad(O);
	}

	double* dist = new double [n_step];
	double* dist_err = new double [n_step];

	double sum = 0;
	double sum2 = 0;

	for(int i=0; i<n_step; i++){
		sum = 0;
		sum2 = 0;
		for(int j=0; j<n_cell; j++){
			sum += dist_quad[i*n_cell + j];
			sum2 += pow(dist_quad[i*n_cell + j],2);
		}
		sum/=n_cell;			//distanza quadratica media
		sum2/=n_cell;			//media del quadrato della distanza quadratica
		dist[i] = sqrt(sum);
		if(i!=0)
			dist_err[i] = sqrt((sum2-pow(sum,2))/i);
		dist_err[i]*=0.5*(1/sqrt(sum));				//propagazione degli errori: incertezza sulla radice del valor medio
	}

	dist_err[0] = 0;

	Stampa("dist_lattice.txt", dist, n_step);
	Stampa("dist_err_lattice.txt", dist_err, n_step);

	//******************************************************//
	
	//---Esempio di RW nel reticolo---//

	int dice = int(rnd.Rannyu(0,n_cell));

	cout<<"Sampled RW on a lattice: "<<dice<<endl; 

	double* X = new double [n_step];
	double* Y = new double [n_step];
	double* Z = new double [n_step];

	for(int i=0; i<n_step; i++){
		X[i] = RW[dice][i].GetX();
		Y[i] = RW[dice][i].GetY();
		Z[i] = RW[dice][i].GetZ();
	}

	Stampa("X.txt", X, n_step);
	Stampa("Y.txt", Y, n_step);
	Stampa("Z.txt", Z, n_step);

	//******************************************************//


	//---Spazio continuo---//

	for(int i=0; i<n_cell; i++){			//inizializzio tutte le posizioni a (0,0,0)
		for(int j=0; j<n_step; j++)
			RW[i][j] = O;
	}

	double theta = 0;					//variabili di appoggio
	double phi = 0;
	sum_x = 0;
	sum_y = 0;
	sum_z = 0;

	for(int i=0; i<n_cell; i++){				//ciclo sulle n_cell RW
		sum_x = 0;
		sum_y = 0;
		sum_z = 0;
		for(int j=1; j<n_step; j++){			//ciclo sugli n_step di ogni RW
			theta = acos(1-2*rnd.Rannyu());		//theta distribuita come (1/2)*sin(theta)
			phi = rnd.Rannyu(0,2*M_PI);		//phi distribuita uniformemente su [0,2pi)
			sum_x += sin(theta)*cos(phi);
			sum_y += sin(theta)*sin(phi);
			sum_z += cos(theta);
			RW[i][j].SetX(sum_x);
			RW[i][j].SetY(sum_y);
			RW[i][j].SetZ(sum_z);
		}
	}

	for(int l=0; l++; l<N)		//vettore delle distanze quadratiche dall'origine
		dist_quad[l] = 0;

	for(int i=0; i<n_step; i++){					//Distanza allo step i-esimo...
		for(int j=0; j<n_cell; j++)				//...nella j-esima RW
			dist_quad[i*n_cell + j] = RW[j][i].Distanza_quad(O);
	}

	for(int l=0; l++; l<n_step){		//inizializzo i vettori usati in precedenza
		dist[l] = 0;
		dist_err[l] = 0;
	}

	sum = 0;
	sum2 = 0;

	for(int i=0; i<n_step; i++){
		sum = 0;
		sum2 = 0;
		for(int j=0; j<n_cell; j++){
			sum += dist_quad[i*n_cell + j];
			sum2 += pow(dist_quad[i*n_cell + j],2);
		}
		sum/=n_cell;			//distanza quadratica media
		sum2/=n_cell;			//media del quadrato della distanza quadratica
		dist[i] = sqrt(sum);
		if(i!=0)
			dist_err[i] = sqrt((sum2-pow(sum,2))/i);
		dist_err[i]*=0.5*(1/sqrt(sum));				//propagazione degli errori: incertezza sulla radice del valor medio
	}

	dist_err[0] = 0;

	Stampa("dist_continuum.txt", dist, n_step);
	Stampa("dist_err_continuum.txt", dist_err, n_step);

	//******************************************************//
	
	//---Esempio di RW nel continuo---//

	dice = int(rnd.Rannyu(0,n_cell));

	cout<<"Sampled RW in continuum: "<<dice<<endl;

	for(int l=0; l<n_step; l++){
		X[l] = 0;
		Y[l] = 0;
		Z[l] = 0;
	}

	for(int i=0; i<n_step; i++){
		X[i] = RW[dice][i].GetX();
		Y[i] = RW[dice][i].GetY();
		Z[i] = RW[dice][i].GetZ();
	}

	Stampa("Xc.txt", X, n_step);
	Stampa("Yc.txt", Y, n_step);
	Stampa("Zc.txt", Z, n_step);

	//******************************************************//

	for(int l=0; l<n_cell; l++)
		delete[] RW[l];
	
	delete[] RW;

	delete[] dist_quad;
	delete[] dist;
	delete[] dist_err;

	delete[] X;
	delete[] Y;
	delete[] Z;

	rnd.SaveSeed();		//Salvo il seme per riproducibilità

	return 0;
}
