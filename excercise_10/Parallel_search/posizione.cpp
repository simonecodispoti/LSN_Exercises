#include "posizione.h"

Posizione::Posizione(){
	m_x = 0.;
	m_y = 0.;
	m_z = 0.;
}

Posizione::Posizione(double x, double y){
	m_x = x;
	m_y = y;
	m_z = 0.;
}

Posizione::Posizione(double x, double y, double z){
	m_x=x;
	m_y=y;
	m_z=z;
}

double Posizione::Get_X() const{
	return m_x;
}

double Posizione::Get_Y() const{
	return m_y;
}

double Posizione::Get_Z() const{
	return m_z;
}

double Posizione::Get_R() const{
	return sqrt( pow(Get_X(),2) + pow(Get_Y(),2) + pow(Get_Z(),2) );
}

double Posizione::Get_Phi() const{
	return atan(Get_Y()/Get_X());
}

double Posizione::Get_Theta() const{
	return acos(Get_Z()/Get_R());
}

double Posizione::Get_Rho() const{
	return sqrt( pow(Get_X(),2) + pow(Get_Y(),2) );
}

double Posizione::Norm_R2(const Posizione& p) const{
	return sqrt( pow(Get_X()-p.Get_X(),2) + pow(Get_Y()-p.Get_Y(),2) );
}

double Posizione::Norm_Quad_R2(const Posizione& p) const{
	return ( pow(Get_X()-p.Get_X(),2) + pow(Get_Y()-p.Get_Y(),2) );
}

double Posizione::Norm_R3(const Posizione& p) const{
	return sqrt( pow(Get_X()-p.Get_X(),2) + pow(Get_Y()-p.Get_Y(),2) + pow(Get_Z()-p.Get_Z(),2) );
}

double Posizione::Norm_Quad_R3(const Posizione& p) const{
	return ( pow(Get_X()-p.Get_X(),2) + pow(Get_Y()-p.Get_Y(),2) + pow(Get_Z()-p.Get_Z(),2) );
}
