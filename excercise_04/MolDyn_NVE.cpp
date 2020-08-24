/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();            //Inizialization
  int nconf = 1;
  for(int iblock=1; iblock <= nblock; ++iblock){
    Reset(iblock);    //Reset block computations
    for(int istep=1; istep <= nstep; ++istep){
      Move();                  //Move particles with Verlet algorithm
      if(istep%10 == 0){
        Measure();              //Properties measurement
        Accumulate();           //Update block computations
        //ConfXYZ(nconf);       //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblock);   //Output for current block
  }
  ConfFinal();          //Write final configuration to restart ---> config.final
  ConfOld();            //Write the previous time step configuration ---> config.old

  return 0;
}


void Input(void){                   //Prepare all stuff for the simulation

  ifstream ReadInput, ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat");

  ReadInput >> old;

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblock;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl << endl;
  cout << "The program uses data-blocking method to estimate averages and uncertainties " << endl;
  cout << "Numer of blocks = " << nblock << endl;
  cout << "Number of steps in each block = " << nstep << endl << endl;

  ReadInput.close();

//Prepare arrays for measurements
  n_props = 4;  // Number of observables
  it = 0;       // Temperature
  iv = 1;       // Potential energy
  ik = 2;       // Kinetic energy
  ie = 3;       // Total energy

//measurement of g(r)
  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;


if(old==0){

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
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
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;	//velocità del cm con massa = 1
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;	//velocità senza cm

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;		              //sto moltiplicando tutte le velocità per il fattore di scala: rad(3*T)
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);     //posizione al tempo t - dt
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
}

if(old==1){

  //Read old configuration
  cout << "Read old configuration from file config.final " << endl << endl;
  ReadConf.open("config.final");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  //Read previous configuration
  cout << "Read previous configuration from file config.old " << endl << endl;
  ReadConf.open("config.old");
  for (int i=0; i<npart; ++i){
    ReadConf >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadConf.close();

  //Rescaling algorithm

  double vx_m[m_part], vy_m[m_part], vz_m[m_part];
  double sumv2=0;
  Move();     // ---> next configuration and velocities

  for(int i=0; i<npart; i++){               //average velocity estimation
    vx_m[i] = Pbc(x[i] - xold[i])/delta;
    vy_m[i] = Pbc(y[i] - yold[i])/delta;
    vz_m[i] = Pbc(z[i] - zold[i])/delta;
  }

  for (int i=0; i<npart; ++i)              //actual kinetic energy
    sumv2 += (vx_m[i]*vx_m[i] + vy_m[i]*vy_m[i] + vz_m[i]*vz_m[i]);

  double fs = sqrt((3*npart*temp)/sumv2);       //veolcity scale factor
  for (int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

    xold[i] = Pbc(x[i] - vx[i] * delta);	//posizione al tempo t - dt
    yold[i] = Pbc(y[i] - vy[i] * delta);
    zold[i] = Pbc(z[i] - vz[i] * delta);
  }
}

  return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

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


double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}


void Measure(){ //Properties measurement
  double bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0;    //reset observables
  t = 0.0;
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;    //reset the hystogram of g(r)

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

      dx = Pbc( x[i] - x[j] );
      dy = Pbc( y[i] - y[j] );
      dz = Pbc( z[i] - z[j] );

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

//update of the histogram of g(r)
      for(int i=0; i<nbins; i++){
        bin = i*bin_size;
        if(dr >= bin and dr < (bin + bin_size)){
          walker[igofr + i] += 2;
          break;
        }
      }

     if(dr < rcut){
      //Potential energy
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart;                 //Potential energy per particle
    stima_kin = t/(double)npart;                 //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart;  //Temperature
    stima_etot = (t+v)/(double)npart;            //Total energy per particle

    walker[iv] = stima_pot;                      //storage for averages (Data-Blocking)
    walker[ik] = stima_kin;
    walker[ie] = stima_etot;
    walker[it] = stima_temp;

    Epot << stima_pot  << endl;                  //istantaneous outputs
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}


void Reset(int iblock){
  if(iblock == 1){
    for(int i=0; i<n_props; ++i){
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for(int i=0; i<n_props; ++i)
    blk_av[i] = 0;

  blk_norm = 0;
}


void Accumulate(void){    //Update block averages
  for(int i=0; i<n_props; ++i)
    blk_av[i] = blk_av[i] + walker[i];

  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){     //Data-blocking technique: printing results for each block
    
   double r, gdir;
   ofstream Epot, Ekin, Etot, Temp, Gofr, Gave;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    
    Epot.open("output_epot.ave",ios::app);
    Ekin.open("output_ekin.ave",ios::app);
    Etot.open("output_etot.ave",ios::app);
    Temp.open("output_temp.ave",ios::app);
    Gofr.open("output.gofr.dat",ios::app);
    
    stima_pot = blk_av[iv]/blk_norm;               //Potential energy per particle
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot = Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_kin = blk_av[ik]/blk_norm;               //Kinetc energy per particle
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin = Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm;               //Total energy per particle
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot = Error(glob_av[ie],glob_av2[ie],iblk);

    stima_temp = blk_av[it]/blk_norm;               //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp = Error(glob_av[it],glob_av2[it],iblk);


//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Kinetic energy per particle
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
//Total energy per particle
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
//Temperature
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;

//g(r): evaluation in each block

    Gofr << iblk << endl;
    for(int k=igofr; k<igofr+nbins; k++){
      r = (k-igofr)*bin_size;                 // distance (from r=0 to r = (nbins-1)*bin_size = 0.5*L-bin_size)
      gdir = blk_av[k]/blk_norm;
      glob_av[k] += gdir;
      glob_av2[k] += gdir*gdir;
      Gofr << setw(wd) << r << setw(wd) << gdir << endl;
    }

//g(r): final estimation with averages and errors

    if(iblk==nblock){    
      Gave.open("output_gofr.ave",ios::app);
      double deltaV = 0;
      double norm = 0;
      for(int k=igofr; k<igofr+nbins; k++){
        r = (k-igofr)*bin_size;
        deltaV = ((4./3.)*M_PI*(pow(r+bin_size,3)-pow(r,3)));     // sphere volume normalization factor
        norm = rho*npart*deltaV;                                  // normalization of gdir
        glob_av[k] /= norm;                                       
        glob_av2[k] /= pow(norm,2);                               // squared norm factor for squared averages
        err_gdir = Error(glob_av[k],glob_av2[k],nblock);
        Gave << setw(wd) << r << setw(wd) << glob_av[k]/(double)nblock << setw(wd) << err_gdir << endl;
      }
      Gave.close();
    }

    cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Gofr.close();
}


void ConfFinal(void){     //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}


void ConfOld(void){     //Write the configuration before the last one
  ofstream WriteConf;

  cout << "Print old configuration to file config.old " << endl << endl;
  WriteConf.open("config.old");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();

  return;
}


void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();

  return;
}


double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}


double Error(double sum, double sum2, int iblk){
  if( iblk == 1 ) return 0.0;
  else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
