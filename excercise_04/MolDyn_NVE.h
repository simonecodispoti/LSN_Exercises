/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//parameters, observables

const int m_props=4;	//number of observables
int n_props;
int iv, ik, it, ie;
double stima_pot, stima_kin, stima_etot, stima_temp;

// averages

double acc, att;

//configuration

const int m_part=108;		//number of molecules
double x[m_part], y[m_part], z[m_part], xold[m_part], yold[m_part], zold[m_part];
double vx[m_part], vy[m_part], vz[m_part];

// thermodynamical state

int npart;					//again, number of molecules
double energy, temp, vol, rho, box, rcut;	//rcut refers to the radius of the cutoff approximation

// simulation

int nstep, ncell, iprint, imeasure, seed;	//total number of integration steps - number of blocks - output frequency - measure frequency - seed for random-gen
double delta;					//time of integration step
bool should_old, should_equi;			//booleans for the options: reading an old config - use equilibration algorithm

//functions

void Input(void);
void Move(void);
void ConfFinal(void);
void ConfPreFinal(void);	//Possibility to save an old configuration
void ConfXYZ(int);
void Measure(void);		//istanteneous values
void Average(void);		//average values
void Eval_ave_err(double* input, double* average, double* error, int n_step, int n_cell);
double Force(int, int);
double Pbc(double);		//Periodic Boundary Conditions
void Equilibrate(void);		//equilibration algorithm

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
