#include "utilities.h"
#include "random.h"
#include "funzione.h"
#include "posizione.h"
using namespace std;

int main(){

	//  || --- | Orbitali s e p: densità di probabilità campionate tramite M(RT)^2 e stima del raggio medio | --- ||  //

	int n_step = pow(10,6);		//numero di step Monte Carlo
	int n_cell = 100;			//numero di blocchi

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

	// --- Uniform sampling: orbitale 1s --- //

	Posizione O;	//Origo initialis

	vector <Posizione> RW;

	for(int i=0; i<n_step; i++)		//RW di n_step: inizializzo tutte le posizioni a (0,0,0)
		RW.push_back(O);

	double r = 1.2;				//passo dell'algoritmo di metropolis: non va scelto nè troppo grande nè troppo piccolo, in modo da garantire un'accettazione del 50% circa
	int ctr = 0;				//contatore delle accettazioni
	const int n_eq = 100;		//numero di step di equilibrazione

	double theta = 0;		//variabili di appoggio
	double phi = 0;
	double x_old = 0;
	double y_old = 0;
	double z_old = 0;
	double x = 0;
	double y = 0;
	double z = 0;
	double q = 0;

	psi_1_0_0* psi = new psi_1_0_0();		//funzione d'onda 1s

	for(int i=1; i<n_step; i++){
		theta = acos(1-2*rnd.Rannyu());		//theta distribuita come (1/2)*sin(theta)
		phi = rnd.Rannyu(0,2*M_PI);			//phi distribuita uniformemente su [0,2pi)
		x = x_old + r*sin(theta)*cos(phi);
		y = y_old + r*sin(theta)*sin(phi);
		z = z_old + r*cos(theta);
		q = (psi -> Eval(x,y,z))/(psi -> Eval(x_old,y_old,z_old));
		if(q > 1 or q > rnd.Rannyu()){
			ctr++;
			x_old = x;
			y_old = y;
			z_old = z;
		}
		if(i > n_eq){				//equilibration steps
			RW[i].Set_X(x_old);
			RW[i].Set_Y(y_old);
			RW[i].Set_Z(z_old);
		}
	}

	cout << "1s orbital - uniform sampling: " << endl;
	cout << "step lenght = " << r << endl;
	cout << "percentage of accepted steps = " << (double(ctr)/n_step)*100 << "%" << endl << endl;		//controllo sulla percentuale di step accettati	

	//Esempio di campionamento per 3D plot

	vector <double> X;
	vector <double> Y;
	vector <double> Z;

	int sample = pow(10,4);

	for(int i=0; i<sample; i++){
		X.push_back(RW[i].Get_X());
		Y.push_back(RW[i].Get_Y());
		Z.push_back(RW[i].Get_Z());
	}

	Print("X_s.txt", X);
	Print("Y_s.txt", Y);
	Print("Z_s.txt", Z);
	
	//Stima del valore medio del raggio (distanza dall'origine)

	vector <double> R;
	for(int i=0; i<n_step; i++)
		R.push_back(RW[i].Norm_R3(O));

	vector <double> media;
	vector <double> errore;
	MC_Mean_Error(R, media, errore, n_step, n_cell);
	Print("Mean_unif_s.txt", media);
	Print("Error_unif_s.txt", errore);


	//---Uniform sampling: orbitale 2p---//

	for(int l=0; l<n_cell; l++)		//inizializzo tutte le posizioni a (0,0,0)
		RW[l] = O;

	r = 2.7;
	ctr = 0;

	theta = 0;
	phi = 0;
	x_old = 0;
	y_old = 0;
	z_old = 0;
	x = 0;
	y = 0;
	z = 0;
	q = 0;

	psi_2_1_0* psi2 = new psi_2_1_0();		//funzione d'onda 2p_z

	for(int i=1; i<n_step; i++){
		theta = acos(1-2*rnd.Rannyu());
		phi = rnd.Rannyu(0,2*M_PI);
		x = x_old + r*sin(theta)*cos(phi);
		y = y_old + r*sin(theta)*sin(phi);
		z = z_old + r*cos(theta);
		q = (psi2 -> Eval(x,y,z))/(psi2 -> Eval(x_old,y_old,z_old));
		if(q > 1 or q > rnd.Rannyu()){
			ctr++;
			x_old = x;
			y_old = y;
			z_old = z;
		}
		if(i > n_eq){
			RW[i].Set_X(x_old);
			RW[i].Set_Y(y_old);
			RW[i].Set_Z(z_old);
		}
	}

	cout << "2p_z orbital - uniform sampling: " << endl;
	cout << "step lenght = " << r << endl;
	cout << "percentage of accepted steps = " << (double(ctr)/n_step)*100 << "%" << endl << endl;		//controllo sulla percentuale di step accettati	

	//Esempio di campionamento per 3D plot

	for(int i=0; i<sample; i++){
		X[i] = RW[i].Get_X();
		Y[i] = RW[i].Get_Y();
		Z[i] = RW[i].Get_Z();
	}

	Print("X_p.txt", X);
	Print("Y_p.txt", Y);
	Print("Z_p.txt", Z);
	
	//Stima del valore medio del raggio (distanza dall'origine)

	for(int i=0; i<n_step; i++)
		R[i] = RW[i].Norm_R3(O);

	MC_Mean_Error(R, media, errore, n_step, n_cell);
	Print("Mean_unif_p.txt", media);
	Print("Error_unif_p.txt", errore);


	// --- Normal sampling: orbitale 1s --- //

	for(int l=0; l<n_cell; l++)		//inizializzo tutte le posizioni a (0,0,0)
		RW[l] = O;

	r = 0.8;
	ctr = 0;

	theta = 0;
	phi = 0;
	x_old = 0;
	y_old = 0;
	z_old = 0;
	x = 0;
	y = 0;
	z = 0;
	q = 0;

	for(int i=1; i<n_step; i++){
		x = rnd.Gauss(x_old,r);
		y = rnd.Gauss(y_old,r);
		z = rnd.Gauss(z_old,r);
		q = (psi -> Eval(x,y,z))/(psi -> Eval(x_old,y_old,z_old));
		if(q > 1 or q > rnd.Rannyu()){
			ctr++;
			x_old = x;
			y_old = y;
			z_old = z;
		}
		if(i > n_eq){
			RW[i].Set_X(x_old);
			RW[i].Set_Y(y_old);
			RW[i].Set_Z(z_old);
		}
	}

	cout << "1s orbital - multivariate normal sampling: " << endl;
	cout << "step lenght = " << r << endl;
	cout << "percentage of accepted steps = " << (double(ctr)/n_step)*100 << "%" << endl << endl;		//controllo sulla percentuale di step accettati	
	
	//Stima del valore medio del raggio (distanza dall'origine)

	for(int i=0; i<n_step; i++)
		R[i] = RW[i].Norm_R3(O);

	MC_Mean_Error(R, media, errore, n_step, n_cell);
	Print("Mean_norm_s.txt", media);
	Print("Error_norm_s.txt", errore);


	//---Normal sampling: orbitale 2p---//

	for(int l=0; l<n_cell; l++)		//inizializzo tutte le posizioni a (0,0,0)
		RW[l] = O;

	r = 1.7;
	ctr = 0;

	theta = 0;
	phi = 0;
	x_old = 0;
	y_old = 0;
	z_old = 0;
	x = 0;
	y = 0;
	z = 0;
	q = 0;

	for(int i=1; i<n_step; i++){
		x = rnd.Gauss(x_old,r);
		y = rnd.Gauss(y_old,r);
		z = rnd.Gauss(z_old,r);
		q = (psi2 -> Eval(x,y,z))/(psi2 -> Eval(x_old,y_old,z_old));
		if(q > 1 or q > rnd.Rannyu()){
			ctr++;
			x_old = x;
			y_old = y;
			z_old = z;
		}
		if(i > n_eq){
			RW[i].Set_X(x_old);
			RW[i].Set_Y(y_old);
			RW[i].Set_Z(z_old);
		}
	}

	cout << "2p_z orbital - multivariate normal sampling: " << endl;
	cout << "step lenght = " << r << endl;
	cout << "percentage of accepted steps = " << (double(ctr)/n_step)*100 << "%" << endl << endl;		//controllo sulla percentuale di step accettati	
	
	//Stima del valore medio del raggio (distanza dall'origine)

	for(int i=0; i<n_step; i++)
		R[i] = RW[i].Norm_R3(O);

	MC_Mean_Error(R, media, errore, n_step, n_cell);
	Print("Mean_norm_p.txt", media);
	Print("Error_norm_p.txt", errore);

	rnd.SaveSeed();	

	return 0;
}
