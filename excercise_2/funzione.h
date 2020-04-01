#ifndef FUNZIONE_H
#define FUNZIONE_H
#include <cmath>
using namespace std;

class FunzioneBase{

	public:

		virtual double Eval(double x) const =0;
};	


//  ||---|funzioni polinomiali|---||  //

class Parabola : public FunzioneBase{
	
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

class Seno: public FunzioneBase{

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


class Coseno: public FunzioneBase{

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

class Seno_smorzato : public FunzioneBase{

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


class Coseno_smorzato : public FunzioneBase{

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

#endif /* FUNZIONE_H */
