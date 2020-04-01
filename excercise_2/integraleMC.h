#ifndef integraleMC_H
#define integraleMC_H
#include "random.h"
#include "funzione.h"
#include <iostream>
#include <cmath>
using namespace std;

class IntegraleMC{
    
	public:
    
		IntegraleMC(Random rnd) {my_rand = rnd;}
		IntegraleMC() {}
		~IntegraleMC() {}
    
		double IntegraleHitOrMiss(double xmin, double xmax, double fmax, int npunti, FunzioneBase* f);		//metodo Hit or Miss semplice
		double IntegraleMedia(double xmin, double xmax, int npunti, FunzioneBase* f);				//metodo della media con sampling uniforme
		double IntegraleMedia(double xmin, double xmax, int npunti, double* vr, FunzioneBase* f);	//metodo della media con importance sampling specificato dal vettore vr
    
	private:
    
		Random my_rand;
};

#endif /* integraleMC_H */
