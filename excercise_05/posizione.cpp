#include "posizione.h"
#include <cmath>
using namespace std;

Posizione::Posizione(){
	m_x=0.;
	m_y=0.;
	m_z=0.;
}

Posizione::Posizione(double x, double y, double z){
	m_x=x;
	m_y=y;
	m_z=z;
}

double Posizione::GetX() const{
	return m_x;
}

double Posizione::GetY() const{
	return m_y;
}

double Posizione::GetZ() const{
	return m_z;
}

double Posizione::GetR() const{
	return sqrt( pow(GetX(),2) + pow(GetY(),2) + pow(GetZ(),2) );
}

double Posizione::GetPhi() const{
	return atan(GetY()/GetX());
}

double Posizione::GetTheta() const{
	return acos(GetZ()/GetR());
}

double Posizione::GetRho() const{
	return sqrt( pow(GetX(),2) + pow(GetY(),2) );
}

double Posizione::Distanza(const Posizione& p) const{
	return sqrt( pow(GetX()-p.GetX(),2) + pow(GetY()-p.GetY(),2) + pow(GetZ()-p.GetZ(),2) );
}

double Posizione::Distanza_quad(const Posizione& p) const{
	return (pow(GetX()-p.GetX(),2) + pow(GetY()-p.GetY(),2) + pow(GetZ()-p.GetZ(),2));
}
