#ifndef VQMC_H
#define VQMC_H
#include <iostream>
#include <cmath>
#include <vector>
#include "random.h"
using namespace std;

//  || --- | potenziali di interazione 1-dimensionali | --- ||  //

class Potential{

	public:

		virtual double Eval(double x) const = 0;
};

class Double_Well : public Potential{

		public:

		Double_Well() {}
		~Double_Well() {}

		virtual double Eval(double x) const {return pow(x,4) - 2.5*pow(x,2);}
};


//  || --- |  dunzioni d'onda test 1-dimensionali | --- ||  //

class Psi_Trial{

	public:

		virtual double Eval(double x) const = 0;
		virtual double SquareMod(double x) const = 0;
		virtual double D2x(double x) const = 0;
};

class Psi_Double_Gauss : public Psi_Trial{

	public:

		Psi_Double_Gauss() {}
		Psi_Double_Gauss(double mu, double sigma) {	m_mu=mu;
													m_sigma=sigma;}
		~Psi_Double_Gauss() {}

		double Get_mu()const {return m_mu;}
		void Set_mu(double mu) {m_mu=mu;}
		double Get_sigma()const {return m_sigma;}
		void Set_sigma(double sigma) {m_sigma=sigma;}

		virtual double Eval(double x) const {return exp(-pow(x-m_mu,2)/(2*pow(m_sigma,2))) + exp(-pow(x+m_mu,2)/(2*pow(m_sigma,2)));}
		virtual double SquareMod(double x) const {return Eval(x)*Eval(x);}
		virtual double D2x(double x) const {return (1/pow(m_sigma,2))*((pow(x-m_mu,2)/pow(m_sigma,2))-1)*exp(-pow(x-m_mu,2)/(2*pow(m_sigma,2))) 
									+ (1/pow(m_sigma,2))*((pow(x+m_mu,2)/pow(m_sigma,2))-1)*exp(-pow(x+m_mu,2)/(2*pow(m_sigma,2)));}

	private:

		double m_mu, m_sigma;
};


class VQ_MC{

    public:
    
        VQ_MC(Random rnd) {m_rnd = rnd;}
        ~VQ_MC() {}

		vector <double> Metropolis_Uniform_Sampling(Psi_Trial* psi, const double start, const double step_size, const int n_step, const int n_eq);
        double Energy_Expected_Value(vector <double> sampling, Potential* pot, Psi_Trial* psi);

	private:

		Random m_rnd;
};

#endif /* VQMC_H */