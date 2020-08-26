#include "VQMC.h"
#include "utilities.h"
using namespace std;

int main(){

	//  || --- | Quantum Variational Monte Carlo | --- ||  //

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

	// --- Optimazion algorithm: SIMULATED ANNEALING --- //

	VQ_MC solver(rnd);
	Double_Well* pot = new Double_Well();				// fix the potential
	Psi_Double_Gauss* psi = new Psi_Double_Gauss();		// fix the trial wave-function

	// Sampling parameters ------------------------------------------------------------------------------------------------------------

	const int sample_step = pow(10,6);		// MC steps for each sampling of the wave-function probability density
	const int sample_eq_step = 1000;		// Equilibration steps for each sampling
	vector <double> sampling;				// Vector of the sampling

	// --------------------------------------------------------------------------------------------------------------------------------

	// Optimization parameters --------------------------------------------------------------------------------------------------------

	const double mu_in = sqrt(5)/2;			// fix the initial conditions
	const double sigma_in = 0.3;
	double mu_step = 0.25;					// fix the initial Metropolis sampling step
	double sigma_step = 0.25;

	vector <int> T;							// setting a temperature range for the cooling
	for(float i=10; i>=0.5; i-=0.5)
		T.push_back(i);

	int n_step = 10;						// setting the number of Metropolis moves for each new cooled system
	int equi_step = 10;

	// --------------------------------------------------------------------------------------------------------------------------------

	// SIMULATED ANNEALING ALGORITHM: changes only in the above lines !!!

	double mu_old = mu_in;					// support variables
	double sigma_old = sigma_in;
	double mu = 0;
	double sigma = 0;
	double energy_old = 0;
	double energy_new = 0;
	double temp = 0;
	double beta = 0;
	double p = 0;
	int ctr = 0;

	for(int i=0; i<T.size(); i++){
		temp = T[i];
		beta = 1/(double)temp;
		mu_step -= 0.01;
		sigma_step -= 0.01;
		ctr = 0;

		for(int j=0; j<equi_step; j++){
			// Old configuration
			psi -> Set_mu(mu_old);
			psi -> Set_sigma(sigma_old);
			sampling = solver.Metropolis_Uniform_Sampling(psi, mu_old, 5*sigma_old, sample_step, sample_eq_step);
			energy_old = solver.Energy_Expected_Value(sampling, pot, psi);
			// New configuration
			mu = rnd.Rannyu(mu_old - 0.5*mu_step, mu_old + 0.5*mu_step);
			sigma = rnd.Rannyu(sigma_old - 0.5*sigma_step, sigma_old + 0.5*sigma_step);
			psi -> Set_mu(mu);
			psi -> Set_sigma(sigma);
			sampling = solver.Metropolis_Uniform_Sampling(psi, mu, 5*sigma, sample_step, sample_eq_step);
			energy_new = solver.Energy_Expected_Value(sampling, pot, psi);
			// Boltzmann acceptance
			if(energy_old > energy_new){
				mu_old = mu;
				sigma_old = sigma;
			}
			continue;
			p = exp(-beta*(energy_new - energy_old));
			if(p > rnd.Rannyu()){
				mu_old = mu;
				sigma_old = sigma;
			}
		}

		if(i==T.size()) n_step += 10;
		for(int j=0; j<n_step; j++){
			// Old configuration
			psi -> Set_mu(mu_old);
			psi -> Set_sigma(sigma_old);
			sampling = solver.Metropolis_Uniform_Sampling(psi, mu_old, 5*sigma_old, sample_step, sample_eq_step);
			energy_old = solver.Energy_Expected_Value(sampling, pot, psi);
			// New configuration
			mu = rnd.Rannyu(mu_old - 0.5*mu_step, mu_old + 0.5*mu_step);
			sigma = rnd.Rannyu(sigma_old - 0.5*sigma_step, sigma_old + 0.5*sigma_step);
			psi -> Set_mu(mu);
			psi -> Set_sigma(sigma);
			sampling = solver.Metropolis_Uniform_Sampling(psi, mu, 5*sigma, sample_step, sample_eq_step);
			energy_new = solver.Energy_Expected_Value(sampling, pot, psi);
			// Boltzmann acceptance
			if(energy_old > energy_new){
				mu_old = mu;
				sigma_old = sigma;
				ctr ++;
			}
			continue;
			p = exp(-beta*(energy_new - energy_old));
			if(p > rnd.Rannyu()){
				mu_old = mu;
				sigma_old = sigma;
				ctr ++;
			}
		}

		cout << "Block: " << i << endl;
		cout << "Temperature: " << T[i] << endl;
		cout << "Acceptance rate: " << ((double)ctr/n_step)*100 << "%" << endl;
		cout << "Optimal parameters: " << "mu = " << mu_old << "    " << "sigma = " << sigma_old << endl;
		cout << "Optimum: " << energy_old << endl;
		cout << "-------------------------------------------------------" << endl;
	}

	// ------------------------------------------------------------------------------------------------------

	// COMMENT THE FOLLOWING SECTION IF INTERESTED ONLY IN THE SEARCH OF OPTIMAL PARAMETERS !!!

	// --- get the probability density for the best values --- //

	psi -> Set_mu(mu_old);
	psi -> Set_sigma(sigma_old);
	sampling = solver.Metropolis_Uniform_Sampling(psi, mu_old, 8*sigma_old, sample_step, sample_eq_step);
	Print("psi2.hist", sampling);

	// --- get the best estimation with data-blocking --- //

	int n_blocks = 100;
	int n_steps = 1000;
	vector <double> energy;

	for(int i=0; i<n_steps; i++){
		sampling = solver.Metropolis_Uniform_Sampling(psi, mu_old, 8*sigma_old, sample_step, sample_eq_step);
		energy.push_back(solver.Energy_Expected_Value(sampling, pot, psi));
	}

	vector <double> average;
	vector <double> error;
	MC_Mean_Error(energy, average, error, n_steps, n_blocks);
	Print("energy.ave", average);
	Print("energy.error", error);

	rnd.SaveSeed();

	return 0;
}