#include "mpi.h"
#include "TSP.h"
#include "utilities.h"
using namespace std;

//  || --- | TSP - GA solver using parallel programming | --- ||  //

int main(int argc, char* argv[]){

	// MPI INITIALIZATION -----------------
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// ------------------------------------

    //******************************PARALLEL_RANDOM_GEN******************************//
	
	// We initiliaze the random generator depending on the rank of the process !!!

	Random rnd;
	int seed[4];
	int p1[size], p2[size];

	ifstream Primes("Primes");
	if (Primes.is_open()){
		for(int i=0; i<size; i++)
			Primes >> p1[i] >> p2[i];
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
				if( property == "RANDOMSEED" ){
					input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
					for(int i=0; i<size; i++)
						if(i==rank) rnd.SetRandom(seed, p1[i], p2[i]);
				}
		}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	//******************************PARALLEL_RANDOM_GEN******************************//

	// --- Settings --- //

	// THESE ARE THE ONLY PARAMETERS THAT CAN BE CHANGED -------------------------
	const int n_cities = 32;		// Problem complexity
	const int n_ind = 500;			// Number of individuals in each population
	const int n_gen = 500;			// Number of generations 
	const int n_mig = 2;			// Number of migrations
	const int n_step = n_ind/2;		// Number of iterations in each generation
	// ---------------------------------------------------------------------------

	// CHOOSE A SINGLE CONFIGURATION AND COMMENT THE OTHER

	// The first continent generate the configuration:

	double* X_Adam = new double[n_cities];
	double* Y_Adam = new double[n_cities];

	// --- Circumference configuration --- //

	/*if(rank==0){
		
		double r = 1.;		// Unitary length
		double theta = 0;
		double x = 0;
		double y = 0;

		for(int i=0; i<n_cities; i++){		// Random placement on a circumference
			theta = rnd.Rannyu(0, 2*M_PI);
			x = r*cos(theta);
			y = r*sin(theta);
			X_Adam[i] = x;
			Y_Adam[i] = y;
		}
	}*/

	// --- Square configuration --- //

	if(rank==0){

		double L = 1.;		// Unitary length
		double x = 0;
		double y = 0;

		for(int i=0; i<n_cities; i++){		// Random placement iside a square
			x = rnd.Rannyu(-L,L);
			y = rnd.Rannyu(-L,L);
			X_Adam[i] = x;
			Y_Adam[i] = y;
		}
	}

	// Wait for the first continent to generate the initialize the configuration:
	MPI_Barrier(MPI_COMM_WORLD);
	// Pass the configuration to the other continents:
	MPI_Bcast(&X_Adam[0], n_cities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Y_Adam[0], n_cities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// Create the same progenitor in every continent:
	Individual Adamo(n_cities);
	for(int i=0; i<n_cities; i++){
		Posizione P(X_Adam[i],Y_Adam[i]);
		City Babele(P);
		Adamo.Set_Gene(Babele,i);
	}
	// Print the initial configuration:
	if(rank==0) Adamo.Print_DNA("Adamo.conf");
	// Create the first generation starting from "Adamo":
	vector <Individual> old_gen = Generation_0(Adamo, n_ind, rnd);
	Pop_Sorting(old_gen);
	// Wait for every continent to have its first generation:
	MPI_Barrier(MPI_COMM_WORLD);

	// --- Genetic Algorithm with migrations --- //

	// MPI TAGS-STATUS to send and receive
	int tag_x = 1;
	int tag_y = 2;
	MPI_Status stat;
	// -----------------------------------

	vector <double> Square_dist;		
	vector <double> Square_dist_ave;

	for(int i=0; i<n_gen; i++){			// Generation cicle
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
		Pop_Sorting(next_gen);		 // Fitness and sorting of the new generation
		// Wait for every continent:
		MPI_Barrier(MPI_COMM_WORLD);

		if(i%n_mig==0){

		// --- migration algorithm --- //

		// For every continent generate the continents tag of migration and the arrays to store migrants DNA:
			int cont_1, cont_2;
			double* X_1 = new double[n_cities];
			double* Y_1 = new double[n_cities];
			double* X_2 = new double[n_cities];
			double* Y_2 = new double[n_cities];
		// The first continent determines randomly the continents of the migration:
			if(rank==0){
				do{
					cont_1 = int(rnd.Rannyu(0,size));
					cont_2 = int(rnd.Rannyu(0,size));
				}while(cont_1 == cont_2);
			}
		// Wait for the first continent:
			MPI_Barrier(MPI_COMM_WORLD);
		// Send the decision:
			MPI_Bcast(&cont_1, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&cont_2, 1, MPI_INT, 0, MPI_COMM_WORLD);
		// Get the migrants and their genes:
			if(rank==cont_1){
				for(int l=0; l<n_cities; l++){
					X_1[l] = next_gen[0].Get_Gene(l).Get_pos().Get_X(); 
					Y_1[l] = next_gen[0].Get_Gene(l).Get_pos().Get_Y();
				}
			}
			if(rank==cont_2){
				for(int l=0; l<n_cities; l++){
					X_2[l] = next_gen[0].Get_Gene(l).Get_pos().Get_X(); 
					Y_2[l] = next_gen[0].Get_Gene(l).Get_pos().Get_Y();
				}
			}
		// Wait for every continent:
			MPI_Barrier(MPI_COMM_WORLD);
		// Optimal individuals migration ---> cont_1 sends to cont_2 and cont_2 sends to cont_1:
			if(rank==cont_1){
				//double* X_best = new double [n_cities];
				//double* Y_best = new double [n_cities];
            	MPI_Send(&X_1[0], n_cities, MPI_DOUBLE, cont_2, tag_x, MPI_COMM_WORLD);
            	MPI_Recv(&X_2[0], n_cities, MPI_DOUBLE, cont_2, tag_x, MPI_COMM_WORLD, &stat);
				//for(int l=0; l<n_cities; l++) X_best[l] = X_2[l];
            	MPI_Send(&Y_1[0], n_cities, MPI_DOUBLE, cont_2, tag_y, MPI_COMM_WORLD);
            	MPI_Recv(&Y_2[0], n_cities, MPI_DOUBLE, cont_2, tag_y, MPI_COMM_WORLD, &stat);
				//for(int l=0; l<n_cities; l++) Y_best[l] = Y_2[l];
				for(int l=0; l<n_cities; l++){
					Posizione Pos(X_2[l], Y_2[l]);
					City city(Pos);
 					next_gen[0].Set_Gene(city, l);
				}
			}
			if(rank==cont_2){
				//double* X_best = new double [n_cities];
				//double* Y_best = new double [n_cities];
            	MPI_Send(&X_2[0], n_cities, MPI_DOUBLE, cont_1, tag_x, MPI_COMM_WORLD);
            	MPI_Recv(&X_1[0], n_cities, MPI_DOUBLE, cont_1, tag_x, MPI_COMM_WORLD, &stat);
				//for(int l=0; l<n_cities; l++) X_best[l] = X_1[l];
            	MPI_Send(&Y_2[0], n_cities, MPI_DOUBLE, cont_1, tag_y, MPI_COMM_WORLD);
            	MPI_Recv(&Y_1[0], n_cities, MPI_DOUBLE, cont_1, tag_y, MPI_COMM_WORLD, &stat);
				//for(int l=0; l<n_cities; l++) Y_best[l] = Y_1[l];
				for(int l=0; l<n_cities; l++){
					Posizione Pos(X_1[l], Y_1[l]);
					City city(Pos);
 					next_gen[0].Set_Gene(city, l);
				}
			}
		}

		Pop_Sorting(next_gen);
		Square_dist.push_back(next_gen[0].Get_Fitness());	// Best quadratic distance of the generation
		old_gen = next_gen;			// saving the actual generation to restart from the present era
	}

	// Print the best path in every continent ---> MODIFY DEPENDING ON THE NUMBER OF CONTINENTS:
	if(rank==0) old_gen[0].Print_DNA("Best_1.conf");
	if(rank==1) old_gen[0].Print_DNA("Best_2.conf");
	if(rank==2) old_gen[0].Print_DNA("Best_3.conf");
	if(rank==3) old_gen[0].Print_DNA("Best_4.conf");
	// Print the best quadratic distance in every continent ---> MODIFY DEPENDING ON THE NUMBER OF CONTINENTS:
	if(rank==0) Print("Distance_1.best", Square_dist);
	if(rank==1) Print("Distance_2.best", Square_dist);
	if(rank==2) Print("Distance_3.best", Square_dist);
	if(rank==3) Print("Distance_4.best", Square_dist);

    rnd.SaveSeed();
	
	MPI_Finalize();
	
    return 0;
}