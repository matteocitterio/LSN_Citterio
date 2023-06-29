#ifndef __Functions__
#define __Functions__

#include <cmath>
#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;


//forward declaration to avoid circular dependency between header files
class Random;

//Parent class used for every basic function
class Functions {

    public:
        virtual double Evaluate (double x) const =0;        //by defining it virtual, we make sure that every child class needs to implement this method
        virtual double Evaluate (arma::vec v) const = 0;
};

class TrialWaveFunction: public Functions {

    public: 
        TrialWaveFunction(Random *rnd);        //empty constructor
        ~TrialWaveFunction() {;};   //empty destructor

        virtual double Evaluate(arma::vec v) const
        {
            cout << "This WF is for a 1D model, hence it works with doubles not vectors. There is something wrong" << endl;
            return 0;
        };
        virtual double Evaluate(double x) const;
        double EvaluateNoModulus(double x) const;
        double SecondDerivative(double x) const;
        void MetropolisUniform();
        void Equilibrate(int nblocks, int L);
        
        void Set_Mu(double mu) {m_mu = mu;}
        void Set_Sigma (double sigma) {m_sigma = sigma;}
        void Set_Position (double x) {m_position =x;}
        
        double Get_Mu() {return m_mu;}
        double Get_Sigma() {return m_sigma;}
        double Get_Position() {return m_position;}

    private:
        double m_mu, m_sigma, m_position, m_stepLength;
        Random *m_rnd;
};

class DoubleDwellPotential: public Functions {

    public:
        DoubleDwellPotential();        //empty constructor
        ~DoubleDwellPotential() {;};   //empty destructor

        virtual double Evaluate(arma::vec v) const{
            cout << "The potential is 1D, please give a 1D INPUT, not a VECTOR" << endl;
            return 0;
        };
        virtual double Evaluate(double x) const;
        void SetCoeffs(arma::vec coeff);

    private:

        double m_a, m_b;

};

#endif // __Functions__