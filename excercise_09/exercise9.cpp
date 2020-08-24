#include "TSP.h"
#include "utilities.h"
using namespace std;

int main(){

    //******************************RANDOM_GEN******************************//
	
	Random rnd;            // this random generator will be used for GA

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

	Individual Adamo (7);		// Progenitor: 32 cities

	double r = 1.;				// Unitary radius

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
	}

	Adamo.Print_DNA("test.dat");
	/*Adamo.Print_DNA();

	Adamo.Per_Mutation(rnd);
	Adamo.Print_DNA();

	Adamo.Shift_Mutation(rnd);
	Adamo.Print_DNA();

	Adamo.Inversion_Mutation(rnd);
	Adamo.Print_DNA();

	Adamo.Swap_Mutation(rnd);
	Adamo.Print_DNA();*/

	vector <Individual> gen_0 = Generation_0(Adamo, 100, rnd);

    rnd.SaveSeed();

    return 0;
}