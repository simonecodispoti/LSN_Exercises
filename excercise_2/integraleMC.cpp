#include "integraleMC.h"

double IntegraleMC::IntegraleHitOrMiss(double xmin, double xmax, double fmax, int npunti, FunzioneBase* f){

    int Nhit=0;
    
    for(int i=0; i<npunti; i++){
        double x = my_rand.Rannyu(xmin, xmax);
        double y = my_rand.Rannyu(0, fmax);
        
        if(y < f->Eval(x))
            Nhit++;
    }
    
    return (xmax-xmin)*fmax*(double(Nhit)/double(npunti));
}

double IntegraleMC::IntegraleMedia(double xmin, double xmax, int npunti, FunzioneBase* f){
    
    double sum = 0.;
    
    for(int i=0; i<npunti; i++){
        sum += (f->Eval(my_rand.Rannyu(xmin, xmax)));
    }
    
    return ((xmax-xmin)*sum)/npunti;
}

double IntegraleMC::IntegraleMedia(double xmin, double xmax, int npunti, double* vr, FunzioneBase* f){
    
    double sum = 0.;
    
    for(int i=0; i<npunti; i++){
        sum += (f->Eval(vr[i]));
    }
    
    return ((xmax-xmin)*sum)/npunti;
}
