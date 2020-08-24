#ifndef UTILITIES_H
#define UTILITIES_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

/* || Library for data analysis and generic utilitities ||

Note: The type T cannot be completely arbitrary! For functions like " MC_Mean_Error " 
T must be a 'numeric' type of some kynd, e.g. int, float, double ...

*/

//  || --- | Input - Output generic functions | --- ||  //

template <class T> void Print(const char* filename, const vector <T>& v){		// Print to file "filename" a single output

	ofstream out;
	out.open(filename);

	if(out.fail()){
		cerr << "Cannot open: "<< filename << "!" << endl;
		exit(-1);
	}else{
		for(int i=0; i<v.size(); i++)
			out << v[i] << endl;
	}

	out.close();
};

template <class T> void Print(const char* filename, const vector <T>& a, const vector <T>& b){

// Print to file "filename" a double output organised in two columns

	ofstream out;
	out.open(filename);

	if(out.fail()){
		cerr << "Cannot open: "<< filename << "!" << endl;
		exit(-2);
	}else{
		for(int i=0; i<a.size(); i++)
			out << setw(12) << a[i] << setw(12) << b[i] << endl;
	}

	out.close();
};

template <class T> void Print(const char* filename, const vector <T>& a, const vector <T>& b, const vector <T>& c){

// Print to file "filename" a triple output organised in three columns

	ofstream out;
	out.open(filename);

	if(out.fail()){
		cerr << "Cannot open: "<< filename << "!" << endl;
		exit(-3);
	}else{
		for(int i=0; i<a.size(); i++)
			out << setw(12) << a[i] << setw(12) << b[i] << setw(12) << c[i] << endl;
	}

	out.close();
};

template <class T> void Print(const vector <T>& v){          // Print to terminal

    for(int i=0; i<v.size(); i++) cout << v[i] << endl;

};

template <class T> vector <T> Read(const char* filename){        // Read from file "filename"

	vector <T> v;

	ifstream in;
	in.open(filename);

	if(in.fail()){
		cerr << "Cannot open: "<< filename << endl;
		exit(-4);
	}else{
		while(!in.eof()){
			T data = 0;
			in >> data;
			v.push_back(data);
		}	
	}
	v.pop_back();		// delete the zero in the end of vector

	in.close();

	return v;
};


//  || --- | Data Analysis | --- ||  //

template <class T> void MC_Mean_Error(const vector <T>& input, vector <T>& average, vector <T>& error, const int n_step, const int n_cell){

	/*
		Data-Blocking technique for the estimation of average values and errors in MC simulations:

		- input:		the vector of type <T> containing datas from a simulation;
		- average:		the vector of type <T> in wich you want to store averages;
		- error:		the vector of type <T> in wich you want to store errors;
		- n_step:		the number of total steps of the simulation (can be different from the total size of "input")
		- n_cell:		the number of blocks in wich datas are divided;
	*/

	int l = n_step/n_cell;
	
	vector <T> ave;
	vector <T> ave2;
	vector <T> sum2_prog;

	for(int i=0; i<n_cell; i++){
		ave.push_back(0);
		ave2.push_back(0);
		sum2_prog.push_back(0);
		average.push_back(0);
		error.push_back(0);
	}

	for(int i=0; i<n_cell; i++){
		double sum = 0;
			for(int j=0; j<l; j++){
				int pos = j+i*l;
				sum += input[pos];
			}
		ave[i] = sum/l;             // averages for each block
		ave2[i] = pow(ave[i],2);   // squared averages for each block
	}

	for(int i=0; i<n_cell; i++){
		for(int j=0; j<i+1; j++){
			average[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
	average[i] /= (i+1);
	sum2_prog[i] /= (i+1);

	if(i!=0)
		error[i] = sqrt((sum2_prog[i]-pow(average[i],2))/i);		// statistical uncertainties
	}

	error[0] = 0;		// no statistical uncertainty after one block
};

template <class T> double Chi2(const vector <T>& observed, const vector <T>& expected){       // Chi2 test
	
	double Chi2 = 0;
    if(observed.size() != expected.size()){
        cerr << "Chi 2 test failed: The size of observations vector must be equal to the size fo expections vector!" << endl;
        exit(-2);
    }else{
       	for(int i=0; i<expected.size(); i++) Chi2 += (pow(observed[i]-expected[i],2))/expected[i];	

	    return Chi2; 
    }
};

#endif /* UTILITIES_H */