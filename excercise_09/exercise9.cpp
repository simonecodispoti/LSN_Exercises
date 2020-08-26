#include "TSP.h"
#include "utilities.h"
using namespace std;

int main(){

	//  || --- | TSP - GA solver  | --- ||  //

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

	// --- Settings --- //

	// THESE ARE THE ONLY PARAMETERS THAT CAN BE CHANGED
	int n_cities = 32;				// Problem complexity
	Individual Adamo(n_cities);		// Progenitor
	int n_ind = 500;				// Number of individuals in each population
	int n_gen = 500;				// Number of generations 
	int n_step = n_ind/2;			// Number of iterations in each generation
	// ---------------------------------------------------------------------------

	// CHOOSE A SINGLE CONFIGURATION AND COMMENT THE OTHER

	// --- Circumference configuration --- //

	/*double r = 1.;		// Unitary radius

	double theta = 0;
	double x = 0;
	double y = 0;

	for(int i=0; i<Adamo.Get_Complexity(); i++){		// Random placement on a circumference
		theta = rnd.Rannyu(0, 2*M_PI);
		x = r*cos(theta);
		y = r*sin(theta);
		Posizione P(x,y);
		City city(P);
		Adamo.Set_Gene(city,i);
	}*/

	// --- Square configuration --- //

	double L =1;		// Unitary length

	double x = 0;
	double y = 0;

	for(int i=0; i<Adamo.Get_Complexity(); i++){		// Random placement iside a square
		x = rnd.Rannyu(-L,L);
		y = rnd.Rannyu(-L,L);
		Posizione P(x,y);
		City city(P);
		Adamo.Set_Gene(city,i);
	}

	// --- Genetic Algorithm --- //

	Adamo.Print_DNA("Adamo.conf");		// Print initial random path

	vector <double> Square_dist;
	vector <double> Square_dist_ave;

	vector <Individual> old_gen = Generation_0(Adamo, n_ind, rnd);		// First generation
	Pop_Sorting(old_gen);				// Fitness and sorting to prepare the next generation
	for(int i=0; i<n_gen; i++){			// Generation cicle
		cout << "Generation " << i << endl;
		cout << "--------------------" << endl; 
		vector <Individual> next_gen;	// Next generation
		for(int j=0; j<n_step; j++){	// Mutations step in each generation
			Individual selected_1 = Natural_Selection(old_gen, rnd);	// Selecting two individuals from sorted population
			Individual selected_2 = Natural_Selection(old_gen, rnd);
			double alea = rnd.Rannyu();
			if(alea <= 0.7)				// Crossover probability = 70%
				Crossover(selected_1, selected_2, rnd);
			if(alea <= 0.1){			// Mutations probability = 10%
				selected_1.Swap_Mutation(rnd);			
				selected_2.Swap_Mutation(rnd);
			}
			else if(alea > 0.1 and alea <= 0.2){
				selected_1.Per_Mutation(rnd);			
				selected_2.Per_Mutation(rnd);
			}
			else if(alea > 0.2 and alea <=0.3){
				selected_1.Inversion_Mutation(rnd);			
				selected_2.Inversion_Mutation(rnd);
			}
			else if(alea > 0.3 and alea <= 0.4){
				selected_1.Shift_Mutation(rnd);			
				selected_2.Shift_Mutation(rnd);
			}
			next_gen.push_back(selected_1);			// Adding the individuals to the next generation	
			next_gen.push_back(selected_2);
		}
		Pop_Sorting(next_gen);		// Fitness and sorting of the new generation
		Square_dist.push_back(next_gen[0].Get_Fitness());		// Best quadratic distance of the generation
		double sum = 0;
		for(int l=0; l<(n_ind/2); l++)		// Average quadratic distance for the best half of the population
			sum += next_gen[l].Get_Fitness();
		sum /= (n_ind/2);
		Square_dist_ave.push_back(sum);
		old_gen = next_gen;		// saving the actual generation to restart from the present era
	}

	old_gen[0].Print_DNA("Best.conf");			// Print the best path
	Print("Distance.best", Square_dist);		// Print the best quadratic distance
	Print("Distance.ave", Square_dist_ave);		// Print the best average quadratic distance

    rnd.SaveSeed();

    return 0;
}