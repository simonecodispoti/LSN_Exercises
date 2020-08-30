#ifndef POSIZIONE_H
#define POSIZIONE_H
#include <cmath>
using namespace std;

class Posizione{

	public:

		Posizione();
		Posizione(double x, double y);
		Posizione(double x, double y, double z);
		~Posizione() {};
		
		double Get_X() const;
		void Set_X(double x) {m_x = x;};

		double Get_Y() const;
		void Set_Y(double y) {m_y = y;};

		double Get_Z() const;
		void Set_Z(double z) {m_z = z;};

		double Get_R() const;
		double Get_Phi() const;
		double Get_Theta() const;
		double Get_Rho() const;

		double Norm_R2(const Posizione&) const;
		double Norm_Quad_R2(const Posizione&) const;

		double Norm_R3(const Posizione&) const;
		double Norm_Quad_R3(const Posizione&) const;

	protected:

		double m_x, m_y, m_z;
};

#endif /* POSIZIONE_H */
