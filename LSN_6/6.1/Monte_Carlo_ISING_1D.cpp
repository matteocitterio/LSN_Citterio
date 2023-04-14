/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  Thermalize();
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  cout << "Seed: " << seed << endl;

  return 0;
}

void Thermalize (void){                     //The first block is used for the equilibration but no measurements are made
  for (int istep=1; istep<=nstep; ++istep){
    Move(metro);
  }
}


void Input(void)
{
  ifstream ReadInput, ReadConf, ReadVelocity, Seed;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2;
   Primes.close();

   ReadInput.open("input.dat");
   ReadInput >> restart;

   if (restart) Seed.open("seed.out");
   else Seed.open("seed.in");
   Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed, p1, p2);
   Seed.close();

   ReadInput >> temp;
   beta = 1.0 / temp;
   cout << "Temperature = " << temp << endl;

   ReadInput >> nspin;
   cout << "Number of spins = " << nspin << endl;

   ReadInput >> J;
   cout << "Exchange interaction = " << J << endl;

   ReadInput >> h;
   cout << "External field = " << h << endl;

   ReadInput >> metro; // if=1 Metropolis else Gibbs

   ReadInput >> nblk;

   ReadInput >> nstep;

   if (metro == 1)
    cout << "The program perform Metropolis moves" << endl;
   else
    cout << "The program perform Gibbs moves" << endl;
   cout << "Number of blocks = " << nblk << endl;
   cout << "Number of steps in one block = " << nstep << endl;

   ReadInput.close();

   // Prepare arrays for measurements
   iu = 0; // Energy
   ic = 1; // Heat capacity
   im = 2; // Magnetization
   ix = 3; // Magnetic susceptibility

   n_props = 4; // Number of observables

   // initial configuration with restart option
   if (restart){
    ReadConf.open("config.final");
    for (int i = 0; i < nspin; ++i)
      ReadConf >> s[i];
   }
   else{
    for (int i = 0; i < nspin; ++i){
      if (rnd.Rannyu() >= 0.5)
        s[i] = 1;
      else
        s[i] = -1;
    }
   }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i){
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    sm = s[o];

    if(metro==1) {            //Run metropolis algorithm

      //compute ene_old
      energy_old = Boltzmann(sm, o);
      //compute ene_new
      energy_new = Boltzmann(-1*sm, o);
      // compute acceptance = exp(-\beta(ene_new-ene_old))
      // calculate probability of flip min(1, acceptance)
      p = min(1., exp(-1 * beta * (energy_new - energy_old)));
      // calculate probability of flip min(1, acceptance)
      // then we need to apply this probability: sample r\in(0,1), if r<Acceptance we flip the spin, else leave it alone
      if (p > rnd.Rannyu()){
        s[o] = -1 * sm;
        accepted += 1;
      }
      attempted+=1;
    }
    
    else  { // Gibbs sampling

      energy_up = Boltzmann(1, o);
      energy_down = Boltzmann(-1, o);
      p = exp(-beta * (energy_down - energy_up));
      if (rnd.Rannyu() < 1 / (1 + p)){
          s[o] = 1;
          accepted++;
      }
      else{
          s[o] = -1;
          accepted++;
      }
      attempted++;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
     // INCLUDE YOUR CODE HERE
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
  // INCLUDE YOUR CODE HERE
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    if (iblk == 1){
     Ene.open("output.ene.0", ios::out);
     Heat.open("output.heat.0", ios::out);
     Mag.open("output.mag.0", ios::out);
     Chi.open("output.chi.0", ios::out);
    }
    else {
     Ene.open("output.ene.0", ios::app);
     Heat.open("output.heat.0", ios::app);
     Mag.open("output.mag.0", ios::app);
     Chi.open("output.chi.0", ios::app);
    }
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    stima_c = pow(beta,2) * (blk_av[ic] / blk_norm - pow(blk_av[iu] / blk_norm, 2)) / (double) nspin; // Heat capacity
    glob_av[ic] += stima_c;
    glob_av2[ic] += stima_c * stima_c;
    err_c = Error(glob_av[ic], glob_av2[ic], iblk);
    Heat << ' ' << iblk << ' ' << stima_c << ' ' << glob_av[ic] / (double)iblk << ' ' << err_c << endl;
    Heat.close();

    stima_m = blk_av[im] / blk_norm / (double)nspin; // Energy
    glob_av[im] += stima_m;
    glob_av2[im] += stima_m * stima_m;
    err_m = Error(glob_av[im], glob_av2[im], iblk);
    Mag << ' ' << iblk << ' ' << stima_m << ' ' << glob_av[im] / (double)iblk << ' ' << err_m << endl;
    Mag.close();

    stima_x = beta * blk_av[ix] / blk_norm / (double)nspin; // Susceptibility
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x * stima_x;
    err_x = Error(glob_av[ix], glob_av2[ix], iblk);
    Chi << ' ' << iblk << ' ' << stima_x << ' ' << glob_av[ix] / (double)iblk << ' ' << err_x << endl;
    Chi.close();

    // INCLUDE YOUR CODE HERE

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
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