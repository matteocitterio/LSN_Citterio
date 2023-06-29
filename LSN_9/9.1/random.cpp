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
#include <cmath>
#include <cstdlib>
#include <armadillo>
#include "random.h"
#include "functions.h"

using namespace std;
using namespace arma;

Random :: Random(){}
// Default constructor, does not perform any action

Random :: ~Random(){}
// Default destructor, does not perform any action

void Random :: SaveSeed(){
   // This function saves the current state of the random number generator to a file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << "RANDOMSEED" <<" "<< l1 << " " << l2 << " " << l3 << " " << l4 << endl;
      ;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   // This function generates a random number from a Gaussian distribution with given mean and sigma
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Exp(double lambda){
   double y = Rannyu();
   return -(log(1 - y)) / (lambda);
}

double Random :: Lorentz(double mu, double gamma){
   double y = Rannyu();
   return mu + gamma * tan(M_PI * (y - 0.5));
}

double Random :: Theta(){
   //samples an angle [0, pi/2]
   double a = Rannyu();
   double b = Rannyu();
   if (b==0){
      return 0;
   }
   if (sqrt(pow(a,2)+pow(b,2))<1){
      return (atan(a/b));
   }
   else{
      return Theta();
   }
}

double Random :: Rannyu(double min, double max){
   // This function generates a random number in the range [min, max]
   return min+(max-min)*Rannyu();
}

double Random ::d_prob(double a, double b, double f_max, Functions *f){
   double x = Rannyu(a, b);
   double r = Rannyu();
   if (r < ((f->Evaluate(x)) / f_max))
   {
      return x;
   }
   else
   {
      return d_prob(a, b, f_max, f);
   }
}

void Random :: Step(vec &r){

   if (r.size() != 3){
      throw invalid_argument("Input vector must have size 3");
   }

   int dim = (int) Rannyu(0,3);
   int s = (int)Rannyu(0, 2) * 2 - 1;

   r[dim] += s;

}

void Random :: CStep (vec &r){
   if (r.size() != 3){
      throw invalid_argument("Input vector must have size 3");
   }

   double x = Gauss(0,1);
   double y = Gauss(0, 1);
   double z = Gauss(0, 1);

   double magnitude = sqrt(pow(x,2)+pow(y,2)+pow(z,2));

   x /= magnitude;
   y /= magnitude;
   z /= magnitude;

   r[0] += x;
   r[1] += y;
   r[2] += z;

}

double Random :: Rannyu(void){
  // This function generates a random number in the range [0,1)
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random ::RandomRoutine(){

  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open())
  {
      Primes >> p1 >> p2;
  }
  else
      cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open())
  {
      while (!input.eof())
      {
         input >> property;
         if (property == "RANDOMSEED")
         {
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            this->SetRandom(seed, p1, p2);
         }
      }
      input.close();
  }
  else
      cerr << "PROBLEM: Unable to open seed.in" << endl;
}

void Random ::SetRandom(int *s, int p1, int p2){
  // This function sets the seed and parameters of the random number generator
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
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
