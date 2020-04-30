#ifndef FUNZIONE_H
#define FUNZIONE_H
#include <cmath>
using namespace std;

class Funzione_R{

	public:

		virtual double Eval(double x) const =0;
};

class Funzione_R3{

	public:
		
		virtual double Eval(double x, double y, double z) const =0;

};


//  ||---|funzioni lineari|---||  //

class Retta : public Funzione_R{
	
	public:

		Retta(double m, double q);
		~Retta() {}

		virtual double Eval(double x) const {return m_m*x + m_q;}

		double Get_m()const {return m_m;}
		void Set_m(double m) {m_m=m;}
		double Get_q()const {return m_q;}
		void Set_q(double q) {m_q=q;}

	private:

		double m_m, m_q;
};


//  ||---|funzioni polinomiali|---||  //

class Parabola : public Funzione_R{
	
	public:

		Parabola(double a, double b, double c);
		~Parabola() {}

		virtual double Eval(double x) const {return m_a*x*x + m_b*x + m_c;}

		double Get_a()const {return m_a;}
		void Set_a(double a) {m_a=a;}
		double Get_b()const {return m_b;}
		void Set_b(double b) {m_b=b;}
		double Get_c()const {return m_c;}
		void Set_c(double c) {m_c=c;}

	private:

		double m_a, m_b, m_c;
};


//  ||---|funzioni trigonometriche|---||  //

class Seno: public Funzione_R{

	public:
	
		Seno(double A, double B, double phi);
		~Seno() {}

		virtual double Eval(double x) const {return m_A*sin(m_B*x + m_phi);}

		double Get_A()const {return m_A;}
		void Set_A(double A) {m_A=A;}
		double Get_B()const {return m_B;}
		void Set_B(double B) {m_B=B;}
		double Get_phi()const {return m_phi;}
		void Set_phi(double phi) {m_phi=phi;}

	private:
		
		double m_A, m_B, m_phi;
};


class Coseno: public Funzione_R{

	public:
	
		Coseno(double A, double B, double phi);
		~Coseno() {}

		virtual double Eval(double x) const {return m_A*cos(m_B*x + m_phi);}

		double Get_A()const {return m_A;}
		void Set_A(double A) {m_A=A;}
		double Get_B()const {return m_B;}
		void Set_B(double B) {m_B=B;}
		double Get_phi()const {return m_phi;}
		void Set_phi(double phi) {m_phi=phi;}

	private:
		
		double m_A, m_B, m_phi;
};


//  ||---|funzioni trigonometriche smorzate|---||  //

class Seno_smorzato : public Funzione_R{

	public:

		Seno_smorzato(double A, double B, double C, double phi);
		~Seno_smorzato() {}

		virtual double Eval(double x) const {return m_A*exp(m_B*x)*sin(m_C*x + m_phi);}

		double Get_A()const {return m_A;}
		void Set_A(double A) {m_A=A;}
		double Get_B()const {return m_B;}
		void Set_B(double B) {m_B=B;}
		double Get_C()const {return m_C;}
		void Set_C(double C) {m_C=C;}
		double Get_phi()const {return m_phi;}
		void Set_phi(double phi) {m_phi=phi;}	

	private:
		
		double m_A, m_B, m_C, m_phi;
};


class Coseno_smorzato : public Funzione_R{

	public:

		Coseno_smorzato(double A, double B, double C, double phi);
		~Coseno_smorzato() {}

		virtual double Eval(double x) const {return m_A*exp(m_B*x)*cos(m_C*x + m_phi);}

		double Get_A()const {return m_A;}
		void Set_A(double A) {m_A=A;}
		double Get_B()const {return m_B;}
		void Set_B(double B) {m_B=B;}
		double Get_C()const {return m_C;}
		void Set_C(double C) {m_C=C;}
		double Get_phi()const {return m_phi;}
		void Set_phi(double phi) {m_phi=phi;}	

	private:
		
		double m_A, m_B, m_C, m_phi;
};


//  ||---|Densità di probabilità orbitali idrogenoidi: notazione {n,l,m}|---||  //

// Si usano coordinate cartesiane e unità ridotte (in raggi di Bohr) //

class psi_1_0_0 : public Funzione_R3{

	public:

		psi_1_0_0();
		~psi_1_0_0() {}

		virtual double Eval(double x, double y, double z) const {return exp(-2*sqrt(pow(x,2)+pow(y,2)+pow(z,2)));}
};

class psi_2_1_0 : public Funzione_R3{

	public:

		psi_2_1_0();
		~psi_2_1_0() {}

		virtual double Eval(double x, double y, double z) const {return pow(z,2)*exp(-sqrt(pow(x,2)+pow(y,2)+pow(z,2)));}
};

#endif /* FUNZIONE_H */
