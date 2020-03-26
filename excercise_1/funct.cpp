#include "funct.h"

void Stampa(const char* filename, const double* v, const int dim){

	ofstream out;
	out.open(filename);

	if(out.fail()){
		cerr<<"Cannot open: "<<filename<<endl;
		exit(-1);
	}else{
		for(int i=0; i<dim; i++)
			out<<v[i]<<endl;
	}

	out.close();
}

void MC_MeanProg(double* media, const int n_step, const int n_cell, const double* data){

	if(n_step % n_cell != 0){
		cerr<<"The simulation must be divided into n_cell blocks: n_step must be a multiple of n_cell!"<<endl;
		exit(-2);
	}

	int l = n_step/n_cell;

	double* ave = new double[n_cell];

	for(int i=0; i<n_cell; i++)
		media[i] =0;

	for(int i=0; i<n_cell; i++){
		double sum = 0;
			for(int j=0; j<l; j++){
				int pos = j+i*l;
				sum += data[pos];
			}
		ave[i] = sum/l;			//media per ogni cella
	}

	for(int i=0; i<n_cell; i++){
		for(int j=0; j<i+1; j++)
			media[i] += ave[j];
	media[i]/=(i+1);
	}

	delete[] ave;
}

void MC_ErrProg(double* errore, const int n_step, const int n_cell, const double* data){

	if(n_step % n_cell != 0){
		cerr<<"The simulation must be divided into n_cell blocks: n_step must be a multiple of n_cell!"<<endl;
		exit(-3);
	}

	int l = n_step/n_cell;

	double* ave = new double [n_cell];
	double* ave2 = new double [n_cell];
	double* sum_prog = new double [n_cell];
	double* sum2_prog = new double [n_cell];

	for(int i=0; i<n_cell; i++){
		sum_prog[i] = 0;
		sum2_prog[i] = 0;
	}

	for(int i=0; i<n_cell; i++){
		double sum = 0;
			for(int j=0; j<l; j++){
				int pos = j+i*l;
				sum += data[pos];
			}
		ave[i] = sum/l;			//media per ogni cella
		ave2[i] = pow(ave[i],2);	//valori medi quadratici per ogni cella
	}

	for(int i=0; i<n_cell; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
	sum_prog[i]/=(i+1);
	sum2_prog[i]/=(i+1);

	if(i!=0)
		errore[i] = sqrt((sum2_prog[i]-pow(sum_prog[i],2))/i);		//incertezza statistica
	}

	errore[0] = 0;		//con un solo blocco non ho incertezza statistica

	delete[] ave;
	delete[] ave2;
	delete[] sum_prog;
	delete[] sum2_prog;
}

double Chi2(const double* obs, const double* exp, const int dim){
	
	double Chi2 = 0;

	for(int i=0; i<dim; i++)
		Chi2 += (pow(obs[i]-exp[i],2))/exp[i];	

	return Chi2;
}









