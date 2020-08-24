/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <cstring>

//parameters, observables
const int m_props=1000;
int n_props, it, iv, ik , ie, igofr, nbins;
double bin_size, sd;
double walker [m_props];

// averages
double blk_av[m_props], blk_norm;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, err_pot, stima_kin, err_kin, stima_etot, err_etot, stima_temp, err_temp, err_gdir;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy, temp, vol, rho, box, rcut;

// simulation
int nstep, nblock, seed;
double delta;
bool old;         //boolean to start from old configurations

//functions
void Input(void);
void Move(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void ConfFinal(void);
void ConfOld(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double,double,int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
