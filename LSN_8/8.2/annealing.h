#ifndef __Annealing__
#define __Annealing__

#include <cmath>
#include <limits>
#include "random.h"
#include "lib.h"
#include "functions.h"
#include <armadillo>

using namespace std;
using namespace arma;

class Annealing{

    public:
        Annealing(TrialWaveFunction *TWF, Random *rnd, DoubleDwellPotential *DDP); // empty constructor
        ~Annealing(){;};    //empty destructor

        void MetropolisAnnealing();
        void SetBeta(double Beta) {m_beta = Beta;}
        void SetEnergy(double Energy) { m_energy = Energy; }
        void SetN(int N) {m_N=N;}
        void SetL(int L) {m_L=L;}

        double GetBeta() { return m_beta; }
        double GetEnergy() {return m_energy;}    
        double GetError() {return m_currentError;}
        double GetMu() {return m_f->Get_Mu();}
        double GetSigma() {return m_f->Get_Sigma();} 

    private : 
        double m_beta;
        double m_energy;
        double m_currentEnergy, m_currentError;
        int m_N;
        int m_L;
        TrialWaveFunction *m_f;
        Random *m_rnd;
        DoubleDwellPotential *m_d;
};


#endif // __Annealing__