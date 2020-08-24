#include "VQMC.h"

vector <double> VQ_MC :: Metropolis_Uniform_Sampling(Psi_Trial* psi, const double start, const double step_size, const int n_step, const int n_eq){

	/*
		Sampling of the distribution given by |psi(x)|^2 using M(RT)^2 uniform distribution sampling, with:

		- start:		the starting point used in equilibration;
		- step_size:	the step length of each Metropolis move;
		- n_step: 		the total number of moves (not counting equilibration);
		- n_eq: 		the number of equilibration moves;
	*/

    vector <double> sample;

	for(int i=0; i<n_step; i++)
		sample.push_back(0);

	int ctr = 0;			//counter of acceptance

	double x_old = start;	//Setting starting point
	double x = 0;
	double q = 0;

	for(int i=0; i<n_eq; i++){			//equilibration
		x = m_rnd.Rannyu(x_old - 0.5*step_size, x_old + 0.5*step_size);
		q = (psi -> SquareMod(x))/(psi -> SquareMod(x_old));
		if(q > 1 or q > m_rnd.Rannyu()) x_old = x;
	}

	for(int i=1; i<n_step; i++){
		x = m_rnd.Rannyu(x_old - 0.5*step_size, x_old + 0.5*step_size);
		q = (psi -> SquareMod(x))/(psi -> SquareMod(x_old));
		if(q > 1 or q > m_rnd.Rannyu()){
			ctr++;
			x_old = x;
		}
		sample[i] = x_old;
	}

	cout << "Sampling finished: " << endl;
	cout << "Rate of accepted moves = " << (double(ctr)/n_step)*100 << "%" << endl;
	cout << "---------------------------------------" << endl << endl;

    return sample;
}

double VQ_MC :: Energy_Expected_Value(vector <double> sampling, Potential* pot, Psi_Trial* psi){

	/*
		Evaluation of the local energy of the form: 

			E(x) = (-0.5/Psi_t(x))*(d^2/dx^2(Psi_t(x))) + V(x)

		and estimation of the expected value of H, using MC importance sampling specified by the vector "sampling"
	*/
    
    double sum = 0.;
    int N = sampling.size();
    
    for(int i=0; i<N; i++){
        sum += ( -0.5*(( psi -> D2x(sampling[i]) )/( psi -> Eval(sampling[i]) )) + pot -> Eval(sampling[i]) );
    }
    
    return sum/N;
}