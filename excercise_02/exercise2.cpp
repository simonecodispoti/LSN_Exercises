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

	vector <double> esiti;

	for(int i=0; i<n_cell ;i++)
		esiti.push_back(integratore.IntegraleMedia(0, 1, n_step, Cos));
	
	Print("Integrali.txt", esiti);

	vector <double> media;
	vector <double> errore;
	MC_Mean_Error(esiti, media, errore, n_cell, n_cell);
	
	Print("I_unif.txt", media);
	Print("I_error_unif.txt", errore);

	delete Cos;

	//---Metodo della media: Importance sampling esponenziale---//

	Coseno_smorzato* f = new Coseno_smorzato(1-exp(-M_PI/2), 0.5*M_PI, 0.5*M_PI, 0);
	IntegraleMC integratore2;

	vector <double> random;
	for(int l=0; l<n_step; l++)
			random.push_back(0);

	vector <double> esiti2;

	for(int i=0; i<n_cell ;i++){
		for(int l=0; l<n_step; l++)
			random[l] = 0;
		for(int j=0; j<n_step; j++)
			random[j] = -(2/M_PI)*log(1-(1-exp(-M_PI/2))*rnd.Rannyu());		// Inverse cumulative method for that distribution
		esiti2.push_back(integratore2.IntegraleMedia(n_step, random, f));
	}
	
	Print("Integrali_exp.txt", esiti2);

	vector <double> media2;
	vector <double> errore2;
	MC_Mean_Error(esiti2, media2, errore2, n_cell, n_cell);

	Print("I_exp.txt", media2);
	Print("I_error_exp.txt", errore2);

	delete f;


	//  ||---|Random walks in 3D|---||  //

	//---Reticolo di spaziatura unitaria---//

	double a = 1;		//spaziatura del reticolo
	n_step = 1000;		//numero di passi del camminatore
	n_cell = 10000;		//Numero di RW

	Posizione O;		//origine (0,0,0)

	vector <vector <Posizione> > RW;
	vector <Posizione> init;

	for(int l=0; l<n_step; l++)
		init.push_back(O);

	for(int i=0; i<n_cell; i++)
		RW.push_back(init);

	double p = 0;					//variabili di appoggio
	double sum_x = 0;
	double sum_y = 0;
	double sum_z = 0;

	for(int i=0; i<n_cell; i++){				//ciclo sulle n_cell RW
		sum_x = 0;
		sum_y = 0;
		sum_z = 0;
		for(int j=1; j<n_step; j++){			//ciclo sugli n_step di ogni RW
			p = rnd.Rannyu();
			if(p < 0.5){						//campionamento di variabile aleatoria discreta tramite partizione dell'intervallo [0,1)
				if(p < double(1)/3){
					if(p < double(1)/6){
						sum_x += a;
						RW[i][j].Set_X(sum_x);
						RW[i][j].Set_Y(sum_y);
						RW[i][j].Set_Z(sum_z);
						continue;
					}
				sum_x -= a;
				RW[i][j].Set_X(sum_x);
				RW[i][j].Set_Y(sum_y);
				RW[i][j].Set_Z(sum_z);
				continue;
				}
			sum_y += a;
			RW[i][j].Set_X(sum_x);
			RW[i][j].Set_Y(sum_y);
			RW[i][j].Set_Z(sum_z);
			}
			else if (p > 0.5){
				if(p > double(2)/3){
					if(p > double(5)/6){
						sum_z -= a;
						RW[i][j].Set_X(sum_x);
						RW[i][j].Set_Y(sum_y);
						RW[i][j].Set_Z(sum_z);
						continue;
					}
				sum_z += a;
				RW[i][j].Set_X(sum_x);
				RW[i][j].Set_Y(sum_y);
				RW[i][j].Set_Z(sum_z);
				continue;
				}
			sum_y -= a;
			RW[i][j].Set_X(sum_x);
			RW[i][j].Set_Y(sum_y);
			RW[i][j].Set_Z(sum_z);
			}
		}
	}

	int N = n_step*n_cell;
	vector <double> dist_quad;		//vettore delle distanze quadratiche dall'origine
	for(int l=0; l<N; l++)
		dist_quad.push_back(0);

	for(int i=0; i<n_step; i++){				//Distanza allo step i-esimo...
		for(int j=0; j<n_cell; j++)				//...nella j-esima RW
			dist_quad[i*n_cell + j] = RW[j][i].Norm_Quad_R2(O);
	}

	vector <double> dist;
	vector <double> dist_err;
	for(int l=0; l<n_step; l++){
		dist.push_back(0);
		dist_err.push_back(0);
	}

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

	Print("dist_lattice.txt", dist);
	Print("dist_err_lattice.txt", dist_err);

	//******************************************************//
	
	//---Esempio di RW nel reticolo---//

	int dice = int(rnd.Rannyu(0,n_cell));

	cout << "Sampled RW on a lattice: "<< dice << endl; 

	vector <double> X;
	vector <double> Y;
	vector <double> Z;

	for(int i=0; i<n_step; i++){
		X.push_back(RW[dice][i].Get_X());
		Y.push_back(RW[dice][i].Get_Y());
		Z.push_back(RW[dice][i].Get_Z());
	}

	Print("X.txt", X);
	Print("Y.txt", Y);
	Print("Z.txt", Z);

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
			phi = rnd.Rannyu(0,2*M_PI);			//phi distribuita uniformemente su [0,2pi)
			sum_x += sin(theta)*cos(phi);
			sum_y += sin(theta)*sin(phi);
			sum_z += cos(theta);
			RW[i][j].Set_X(sum_x);
			RW[i][j].Set_Y(sum_y);
			RW[i][j].Set_Z(sum_z);
		}
	}

	for(int l=0; l++; l<N)		//vettore delle distanze quadratiche dall'origine
		dist_quad[l] = 0;

	for(int i=0; i<n_step; i++){			//Distanza allo step i-esimo...
		for(int j=0; j<n_cell; j++)			//...nella j-esima RW
			dist_quad[i*n_cell + j] = RW[j][i].Norm_Quad_R2(O);
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

	Print("dist_continuum.txt", dist);
	Print("dist_err_continuum.txt", dist_err);

	//******************************************************//
	
	//---Esempio di RW nel continuo---//

	dice = int(rnd.Rannyu(0,n_cell));

	cout << "Sampled RW in continuum: "<< dice << endl;

	for(int i=0; i<n_step; i++){
		X[i] = RW[dice][i].Get_X();
		Y[i] = RW[dice][i].Get_Y();
		Z[i] = RW[dice][i].Get_Z();
	}

	Print("Xc.txt", X);
	Print("Yc.txt", Y);
	Print("Zc.txt", Z);

	//******************************************************//

	rnd.SaveSeed();		//Salvo il seme per riproducibilità

	return 0;
}
