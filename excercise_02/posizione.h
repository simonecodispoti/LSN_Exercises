#ifndef POSIZIONE_H
#define POSIZIONE_H
using namespace std;

class Posizione{

	public:

		Posizione();
		Posizione(double x, double y, double z);
		~Posizione() {};
		
		double GetX() const;
		void SetX(double x) {m_x = x;};

		double GetY() const;
		void SetY(double y) {m_y = y;};

		double GetZ() const;
		void SetZ(double z) {m_z = z;};

		double GetR() const;
		double GetPhi() const;
		double GetTheta() const;
		double GetRho() const;

		double Distanza(const Posizione&) const;
		double Distanza_quad(const Posizione&) const;

	protected:

		double m_x, m_y, m_z;
};

#endif /* POSIZIONE_H */
