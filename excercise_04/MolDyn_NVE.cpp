/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>	//srand, rand: to generate random number
#include <iostream>
#include <fstream>
#include <cmath>
#include "MolDyn_NVE.h"

using namespace std;

int main(){

  Input();										//Inizialization
  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
     Move();										//Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%imeasure == 0){
        Measure();									//Properties measurement
        //ConfXYZ(nconf);								//Write actual configuration in XYZ format - Commented to avoid "filesystem full"! 
        nconf += 1;
     }
  }
  cout << endl;
  Average();										//compute average values and errors using blocking method
  ConfFinal();										//Write final configuration to restart
  ConfPreFinal();									//Write olso the old configuration (possibility to rester a simulation)

  return 0;
}


void Input(void){									//Prepare all stuff for the simulation

  ifstream ReadInput, ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential: v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl << endl;

  seed = 1;		//Set seed for random numbers
  srand(seed);		//Initialize random number generator

  //------------------------------------------------------------------------------------//
  
  ReadInput.open("input.dat");

  ReadInput >> should_old;
  if(should_old==true)
    cout<< "Starting from an old configuration config.old " << endl << endl;

  ReadInput >> should_equi;
  if(should_equi==true)
    cout<< "Equilibration of the ensemble option enabled " << endl << endl;

  cout << "Ensemble properties: " << endl << endl;

  ReadInput >> temp;
  cout << "Temperature of the ensemble = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> ncell;
  ReadInput >> imeasure;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations using Verlet algorithm " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  cout << "Number of blocks = " << ncell << endl;
  cout << "Measurement frequency = " << imeasure << endl << endl;

  ReadInput.close();

  //------------------------------------------------------------------------------------//

  //Prepare array for measurements

  iv = 0;	//Potential energy
  ik = 1;	//Kinetic energy
  ie = 2;	//Total energy
  it = 3;	//Temperature
  n_props = 4;	//Number of observables

  //Read initial configuration

  cout << "Reading initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

   //Prepare initial velocities

   //case 1: no old configurations

   if(should_old==false){

     cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
     double sumv[3] = {0.0, 0.0, 0.0};
     for (int i=0; i<npart; ++i){
       vx[i] = rand()/double(RAND_MAX) - 0.5;
       vy[i] = rand()/double(RAND_MAX) - 0.5;
       vz[i] = rand()/double(RAND_MAX) - 0.5;

       sumv[0] += vx[i];
       sumv[1] += vy[i];
       sumv[2] += vz[i];
     }

     for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;	//center mass veolocity for identical particles
     double sumv2 = 0.0, fs;
     for (int i=0; i<npart; ++i){
       vx[i] = vx[i] - sumv[0];						//subtracting center of mass velocity components
       vy[i] = vy[i] - sumv[1];
       vz[i] = vz[i] - sumv[2];

       sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];			
     }
     sumv2 /= (double)npart;

     fs = sqrt(3 * temp / sumv2);		// fs = velocity scale factor 
     for (int i=0; i<npart; ++i){
       vx[i] *= fs;
       vy[i] *= fs;
       vz[i] *= fs;

       xold[i] = Pbc(x[i] - vx[i] * delta);
       yold[i] = Pbc(y[i] - vy[i] * delta);
       zold[i] = Pbc(z[i] - vz[i] * delta);
     }

   }

   //case 2: using an old configuration

   if(should_old==true){

     cout << "Reading old configuration from file config.old " << endl << endl;
     ReadConf.open("config.old");
     for (int i=0; i<npart; ++i){
       ReadConf >> xold[i] >> yold[i] >> zold[i];
       xold[i] = xold[i] * box;
       yold[i] = yold[i] * box;
       zold[i] = zold[i] * box;
     }

   }

   if(should_equi==true)
     Equilibrate();		//Equilibration before starting the true simulation

   return;
}

void Equilibrate(void){

   /*cout << "Equilibration algorithm starting " << endl;
   double prec = 0;
   cout << "Choose a precision to fix thermodynamic phases (note that a too high precision can cause an infinite loop!): " << endl;
   cin >> prec;*/

   double vx_mid[m_part], vy_mid[m_part], vz_mid[m_part];
   double Kin_mid = 0.;
   double T_mid = 0.;
   double scale;

   Move();

   for (int i=0; i<npart; i++){
     vx_mid[i] = Pbc(x[i] - xold[i])/delta;
     vy_mid[i] = Pbc(y[i] - yold[i])/delta;
     vz_mid[i] = Pbc(z[i] - zold[i])/delta;
   }
  
   for (int i=0; i<npart; ++i)
     Kin_mid += 0.5 * (vx_mid[i]*vx_mid[i] + vy_mid[i]*vy_mid[i] + vz_mid[i]*vz_mid[i]);

   Kin_mid/=(double)npart;
   T_mid = (2.0 / 3.0) * Kin_mid;

   scale = sqrt(temp/T_mid);

   for (int i=0; i<npart; i++){
     vx[i] *= scale;
     vy[i] *= scale;
     vz[i] *= scale;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }

  return;
}

void Move(void){ 							//Move particles with Verlet algorithm

  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ 						//Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ 						//Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){						//Compute forces as -Grad_ip V(r)

  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );					// distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8));		// -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){					//Properties measurement

  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0;	//reset observables
  t = 0.0;

  //cycle over pairs of particles

  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

  //Potential energy
       v += vij;
     }
    }          
  }

  //Kinetic energy

  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart;		//Potential energy per particle
    stima_kin = t/(double)npart;		//Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart;		//Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}

void Average(void){

  ifstream Epot, Ekin, Etot, Temp;
  ofstream Epot_ave, Ekin_ave, Etot_ave, Temp_ave, Epot_err, Ekin_err, Etot_err, Temp_err;

  //reading from outputs of Measure()

  Epot.open("output_epot.dat");
  Ekin.open("output_ekin.dat");
  Temp.open("output_temp.dat");
  Etot.open("output_etot.dat");

  int mis = nstep/imeasure;	//number of measurements
 
  double* ep = new double [mis];
  double* et = new double [mis];
  double* ek = new double [mis]; 
  double* T = new double [mis];

  for(int i=0; i<mis; i++){
    Epot >> ep[i];
    Ekin >> ek[i];
    Etot >> et[i];
    Temp >> T[i];
  }

  Epot.close();
  Ekin.close();
  Etot.close();
  Temp.close();

  //using data blocking for each observable

  Epot_ave.open("output_epot_ave.dat");
  Ekin_ave.open("output_ekin_ave.dat");
  Temp_ave.open("output_temp_ave.dat");
  Etot_ave.open("output_etot_ave.dat");
  Epot_err.open("output_epot_err.dat");
  Ekin_err.open("output_ekin_err.dat");
  Temp_err.open("output_temp_err.dat");
  Etot_err.open("output_etot_err.dat");

  double* ep_ave = new double [ncell];
  double* et_ave = new double [ncell];
  double* ek_ave = new double [ncell]; 
  double* T_ave = new double [ncell];
  double* ep_err = new double [ncell];
  double* et_err = new double [ncell];
  double* ek_err = new double [ncell]; 
  double* T_err = new double [ncell];

  Eval_ave_err(ep, ep_ave, ep_err, mis, ncell);
  Eval_ave_err(ek, ek_ave, ek_err, mis, ncell);
  Eval_ave_err(et, et_ave, et_err, mis, ncell);
  Eval_ave_err(T, T_ave, T_err, mis, ncell);

  for(int i=0; i<ncell; i++){
    Epot_ave << ep_ave[i] << endl;
    Ekin_ave << ek_ave[i] << endl;
    Etot_ave << et_ave[i] << endl;
    Temp_ave << T_ave[i] << endl;
    Epot_err << ep_err[i] << endl;
    Ekin_err << ek_err[i] << endl;
    Etot_err << et_err[i] << endl;
    Temp_err << T_err[i] << endl;
  }

  Epot_ave.close();
  Ekin_ave.close();
  Temp_ave.close();
  Etot_ave.close();
  Epot_err.close();
  Ekin_err.close();
  Temp_err.close();
  Etot_err.close();

  delete[] ep;
  delete[] ek;
  delete[] et;
  delete[] T;
  delete[] ep_ave;
  delete[] ek_ave;
  delete[] et_ave;
  delete[] T_ave;
  delete[] ep_err;
  delete[] ek_err;
  delete[] et_err;
  delete[] T_err;

  return;
}

void Eval_ave_err(double* input, double* average, double* error, int n_step, int n_cell){

	int l = n_step/n_cell;
	
	double ave [ncell];
	double ave2 [n_cell];
	double sum2_prog [n_cell];

	for(int i=0; i<n_cell; i++){
		ave[i] = 0;
		ave2[i] = 0;
		sum2_prog[i] = 0;
		average[i] = 0;
		error[i] = 0;
	}

	for(int i=0; i<n_cell; i++){
		double sum = 0;
			for(int j=0; j<l; j++){
				int pos = j+i*l;
				sum += input[pos];
			}
		ave[i] = sum/l;			//averages for each block
		ave2[i] = pow(ave[i],2);	//squared averages for each block
	}

	for(int i=0; i<n_cell; i++){
		for(int j=0; j<i+1; j++){
			average[i] += ave[j];
			sum2_prog[i] += ave2[j];
		}
	average[i]/=(i+1);
	sum2_prog[i]/=(i+1);

	if(i!=0)
		error[i] = sqrt((sum2_prog[i]-pow(average[i],2))/i);		//statistical uncertainty (output)
	}

	error[0] = 0;		//no statistical uncertainty after one block

  return;
}

void ConfFinal(void){			//Write final configuration

  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfPreFinal(void){		//Write the old configuration

  ofstream WriteConf;

  cout << "Print old configuration to file config.old " << endl << endl;
  WriteConf.open("config.old");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){		//Write configuration in .xyz format

  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;

  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){			//Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
