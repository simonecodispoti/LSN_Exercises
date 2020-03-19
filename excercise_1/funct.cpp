#include "funct.h"

void Stampa(const char* filename, const double* v, const int dim){

	ofstream out;
	out.open(filename);

	if(out.fail()){
		cerr<<"Cannot open: "<<filename<<endl;
		exit(-2);
	}else{
		for(int i=0; i<dim; i++)
			out<<v[i]<<endl;
	}
}
