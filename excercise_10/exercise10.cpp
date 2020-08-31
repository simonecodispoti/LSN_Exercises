#include "TSP.h"
#include "utilities.h"
using namespace std;

int main(){

    //******************************RANDOM_GEN******************************//
	
	Random rnd;            // this random generator is going to be used for GA

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

	// --- Settings and annealing schedule --- //

	// THESE ARE THE ONLY PARAMETERS THAT CAN BE CHANGED
	int n_cities = 32;					// Problem complexity
	Individual Initial(n_cities);		// System of cities
	double max_temp = 100;				// Starting temperature
	double min_temp = 0.001;			// Cooling limit
	double scale = max_temp/1000;		// Cooling scale factor
	int n_step = pow(10,4);				// Number of iterations for each temperature
	// -------------------------------------------------------------------------------

	// CHOOSE A SINGLE CONFIGURATION AND COMMENT THE OTHER

	// --- Circumference configuration --- //

	double r = 1.;		// Unitary radius

	double theta = 0;
	double x = 0;
	double y = 0;

	for(int i=0; i<Initial.Get_Complexity(); i++){		// Random placement on a circumference
		theta = rnd.Rannyu(0, 2*M_PI);
		x = r*cos(theta);
		y = r*sin(theta);
		Posizione P(x,y);
		City city(P);
		Initial.Set_Gene(city,i);
	}

	// --- Square configuration --- //

	/*double L = 1.;		// Unitary length

	double x = 0;
	double y = 0;

	for(int i=0; i<Initial.Get_Complexity(); i++){		// Random placement iside a square
		x = rnd.Rannyu(-L,L);
		y = rnd.Rannyu(-L,L);
		Posizione P(x,y);
		City city(P);
		Initial.Set_Gene(city,i);
	}*/

	// --- Simulated Annealing Algorithm --- //

	Initial.Print_DNA("Initial.conf");	// Print initial configuration

	vector <double> Square_dist;		// containers to print informations
	vector <double> Temperature;

	Individual conf_old = Initial;		// Support variables
	Individual conf_new;
	double beta = 0;
	double cost_old = 0;
	double cost_new = 0;
	double p = 0;

	for(double i=max_temp; i>=min_temp; i-=(i*scale)){		// Temperature cicle
		beta = 1/i;
		cout << "Temperature = " << i << endl;
		cout << "----------------------------" << endl; 
		Temperature.push_back(i);
		for(int j=0; j<n_step; j++){			// Mutations step for each cooling
			conf_new = conf_old;				// Saving old configuration for mutation

			double alea = rnd.Rannyu();
			if(alea <= 0.25)					// Mutations probability = 25%
				conf_new.Swap_Mutation(rnd);		
			else if(alea > 0.25 and alea <= 0.50)
				conf_new.Per_Mutation(rnd);
			else if(alea > 0.50 and alea <=0.75)
				conf_new.Inversion_Mutation(rnd);	
			else if(alea > 0.75)
				conf_new.Shift_Mutation(rnd);

			conf_old.Eval_Fitness();			// Fitness cost evaluation
			cost_old = conf_old.Get_Fitness();
			conf_new.Eval_Fitness();
			cost_new = conf_new.Get_Fitness();

			// Boltzmann acceptance
			if(cost_old > cost_new){
				conf_old = conf_new;
				Square_dist.push_back(cost_new);
				continue;
			}
			p = exp(-beta*(cost_new- cost_old));
			if(p > rnd.Rannyu()){
				conf_old = conf_new;
				Square_dist.push_back(cost_new);
			}		

			Square_dist.push_back(cost_old);
		}
	}

	Print("Temperatures.dat", Temperature);		// Print the annealing schedule
	conf_old.Print_DNA("Best.conf");			// Print the best path
	Print("Distance.best", Square_dist);		// Print the best quadratic distance

    rnd.SaveSeed();

    return 0;
}