#include "funzione.h"

Retta::Retta(double m, double q){
	m_m=m;
	m_q=q;
}

Parabola::Parabola(double a, double b, double c){
	m_a=a;
	m_b=b;
	m_c=c;
}

Seno::Seno(double A, double B, double phi){
	m_A=A;
	m_B=B;
	m_phi=phi;
}

Coseno::Coseno(double A, double B, double phi){
	m_A=A;
	m_B=B;
	m_phi=phi;
}

Seno_smorzato::Seno_smorzato(double A, double B, double C, double phi){
	m_A=A;
	m_B=B;
	m_C=C;
	m_phi=phi;
}

Coseno_smorzato::Coseno_smorzato(double A, double B, double C, double phi){
	m_A=A;
	m_B=B;
	m_C=C;
	m_phi=phi;
}

psi_1_0_0::psi_1_0_0() {}

psi_2_1_0::psi_2_1_0() {}
