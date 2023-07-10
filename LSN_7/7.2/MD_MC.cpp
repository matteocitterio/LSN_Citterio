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
#include <vector>
#include "MD_MC.h"

using namespace std;

int main()  {

  Input();                                                                          // Inizialization and input readings  
  Thermalize();                                                                     // Equilibration

  int nconf = 1;

  // Simulation
  for(int iblk=1; iblk <= nblk; iblk++){                                            // Loop over the number of blocks
    Reset(iblk);                                                                    // Reset block averages

    for(int istep=1; istep <= nstep; istep++) {                                     // Loop over the steps per block

      Move();                                                                       // Make the system evolve with the desired algorithm (Metropolis or MD)
      Measure();                                                                    // Measure the relevant quantities
      Accumulate();                                                                 // Update block averages
      if(istep%10 == 0){
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk);                                                                 // Print results for current block & I/O managment
  }

  ConfFinal();                                                                      // Write final configuration
  return 0;
  }

void Thermalize(void){

  /*
  This function takes care of the first thermalization of the system without measuring properties
  */

  cout << "Doing equilibration" << endl;

  for (int istep = 1; istep <= 10000; ++istep) {                                    // The equilibration process takes a block of steps
    Move();
  }
  
}

void Input() {

  /*
  Variables initialization and input readingss
  */

  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  // User friendly messages
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  // Read seed for random numbers
  Primes.open("Primes");
  int p1, p2;
  Primes >> p1 >> p2;
  Primes.close();

  // Read input informations
  ReadInput.open("input.in");

  ReadInput >> phaseName;
  ReadInput >> iNVET;
  ReadInput >> restart;

  if (restart)
    Seed.open("seed.out");
  else
    Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  Seed.close();

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Phase name: " << phaseName << endl;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

  // Tail corrections
  vtail = (8 * M_PI * rho) * ((1 / (9 * pow(rcut, 9))) - (1 / (3 * pow(rcut, 3))));
  ptail = (32 * M_PI * rho) * ((1 / (9 * pow(rcut, 9))) - (1 / (6 * pow(rcut, 3))));

  cout << "U tail correction: " << vtail << endl;
  cout << "P tail correction: " << ptail << endl;

  // Prepare arrays for measurements
  iv = 0;                                                                           // Potential energy
  it = 1;                                                                           // Temperature
  ik = 2;                                                                           // Kinetic energy
  ie = 3;                                                                           // Total energy
  iw = 4;                                                                           // Pressure

  //gr histogram
  bin_size = box /  (2.0 * (double)n_bins);

  // Read initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart)
  {
    ReadConf.open("/config.out");
    ReadVelocity.open("/velocity.out");
    for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
  }
  else 
  {
    ReadConf.open("config.in");
    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i)
    {
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    cout << "velocity scale factor: " << fs << endl << endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();

  for (int i=0; i<npart; ++i)
  {
    if(iNVET)
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
  // Evaluate properties of the initial configuration
  Measure();

  // Print initial values for measured properties
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "Initial potential energy with corrections: " << walker[iv] / (double)npart + vtail << endl;
  cout << "Initial temperature      = " << walker[it] << endl;
  cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;
  cout << "Initial total pressure   = " << walker[iw]<< endl;
  cout << "Initial total pressure with tail corrections  = " << walker[iw] +ptail << endl;

  return;
}

void Move() {

  /*
  Evolves the system according to the algorithm chosen through input.in file (MC | MD)
  */

  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) {                                                                       // Monte Carlo (NVT) move
    for(int i=0; i<npart; ++i) {                                                    // Loop over the number of particles in our simulation box
      o = (int)(rnd.Rannyu()*npart);                                                // Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)

      energy_old = Boltzmann(x[o],y[o],z[o],o);                                     // Compute the energy of the Old configuration

      // Compute the new coordinates of the particle using PBCs
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o],y[o],z[o],o);                                     // Compute the energy of the new configuration
      p = exp(beta*(energy_old-energy_new));                                        // Metropolis acceptance

      if(p >= rnd.Rannyu())  {                                                      // Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted = accepted + 1.0;
      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0;
    }
  } else {                                                                          // Molecular Dynamics (NVE) move
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){                                                     // Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){                                                     // Verlet integration scheme

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

      accepted = accepted + 1.0;
      attempted = attempted + 1.0;
    }
  }
  return;
}

double Boltzmann(double xx, double yy, double zz, int ip) {

  /*
  This function computes the energy contribute to the system of a particle placed in position (xx, yy, zz)
  */

  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip) {                                                                   // Avoid counting self contributes
      // distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Force(int ip, int idir){                                                     // Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );                                                // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8));                       // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Measure(){                                                                     

  /*
  Measurment of relevant quantities
  */

  double v = 0.0, kin=0.0, press = 0.0;
  double vij;
  double pij;
  double dx, dy, dz, dr;
  int bin;

  for (int i = 0; i<n_bins; i++) walker[i+n_props]=0.0;                             // Reset gr histo

  for (int i=0; i<npart-1; ++i){                                                    // Cycle over pairs of particles
    for (int j=i+1; j<npart; ++j) {                                                 // The other particles
      // distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut) {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        v += vij;
        pij = 1.0 / pow(dr, 12) - 0.5 / pow(dr, 6);
        press += pij;       
      }

      if (dr <= (box / 2)) {                                                        // My Gr is built with r \in [0, box/2]
        bin = (int)(n_props + (dr / bin_size));                                     // Choose the correct bin
        walker[bin] += 2;
      }
    }          
  }

  for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  walker[iv] = 4.0 * v;                                                             // Potential energy
  walker[ik] = kin;                                                                 // Kinetic energy
  walker[it] = (2.0 / 3.0) * kin/(double)npart;                                     // Temperature
  walker[ie] = 4.0 * v + kin;                                                       // Total energy;
  walker[iw] = walker[it] * rho + 48 / (3 * vol) * press;                           //pressure

  // Rercord instant values
  instantU.open("./" + phaseName + "/instant.epot", ios::app);
  instantP.open("./" + phaseName + "/instant.pres", ios::app);

  instantU << walker[iv] / (int) npart + vtail << endl;
  instantP << walker[iw] + ptail << endl;

  instantU.close();
  instantP.close();

  return;
}


void Reset(int iblk) {

  /*
  Reset block averages
  */

  if (iblk == 1) {
    for (int i = 0; i < m_props; ++i) {
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for(int i=0; i<m_props; ++i){
     blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}


void Accumulate(void) {

  /*
  Update block averages
  */

  for (int i = 0; i < m_props; ++i){
     blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) {

  /*
  Print results for current block
  */

  ofstream Epot, Press, Gr;

  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted / attempted << endl<< endl;

  if (iblk == 1) {
    Epot.open(phaseName + "/U.dat", ios::trunc);                                   // Rewriting the file for the first block
    Press.open(phaseName + "/P.dat", ios::trunc);
    Gr.open(phaseName + "/Gr.dat", ios::trunc);
  }
  else {
    Epot.open(phaseName + "/U.dat", ios::app);                                      // Appending data for subsequent blocks
    Press.open(phaseName + "/P.dat", ios::app);
    Gr.open(phaseName + "/Gr.dat", ios::app);
  }

  stima_pot = blk_av[iv] / blk_norm / (double)npart + vtail;                        // Potential energy per particle
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot * stima_pot;
  err_pot = Error(glob_av[iv], glob_av2[iv], iblk);

  stima_pres = blk_av[iw] / blk_norm + ptail;                                       // Pressure
  glob_av[iw] += stima_pres;
  glob_av2[iw] += stima_pres * stima_pres;
  err_press = Error(glob_av[iw], glob_av2[iw], iblk);

  // File output
  Epot << " " << iblk << " " << stima_pot << " " << glob_av[iv] / (double)iblk << " " << err_pot << endl;
  Press << " " << iblk << " " << stima_pres << " " << glob_av[iw] / (double)iblk << " " << err_press << endl;

  // Gr histogram for every block
  for (int i = 0; i < n_bins; i++) {                                                // Loop over the bins
    int h = (int)(n_props + i);
    deltaV = 4 * M_PI * (pow((bin_size * (i + 1)), 3) - pow((bin_size * i), 3)) / 3;
    stima_gr = blk_av[h] / blk_norm / (rho * npart * deltaV);
    glob_av[h] += stima_gr;
    glob_av2[h] += stima_gr * stima_gr;
    Gr << " " << iblk << " " << bin_size * i << " " << stima_gr << " " << glob_av[h] / (double)iblk << endl;
  }

  if (iblk == nblk) {                                                               // Print final Gr with errors 
    ofstream finalGr;
    finalGr.open(phaseName + "/finalGr.dat", ios::trunc);
    for (int i = 0; i < n_bins; i++)  {                                             // Loop over the bins
      int h = (int)(n_props + i);
      err_gr[i] = Error(glob_av[h], glob_av2[h], nblk);
      finalGr << " " << bin_size * i << " " << " " << glob_av[h] / (double)iblk << " " << err_gr[i] << endl;
    }
    finalGr.close();
  }

  cout << "----------------------------" << endl << endl;

  Epot.close();
  Press.close();
  Gr.close();
  
}

void ConfFinal()
{
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("config.out");
  WriteVelocity.open("velocity.out");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
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
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if (iblk == 1) return 0.0;
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

