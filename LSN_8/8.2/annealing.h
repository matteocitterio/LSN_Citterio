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

    /*
    This class takes care of the simulated annealing schema for a trial wave function and an input dwell-shaped potential
    */

    public:
        Annealing(TrialWaveFunction *TWF, Random *rnd, DoubleDwellPotential *DDP);  // empty constructor
        ~Annealing(){;};                                                            //empty destructor

        void MetropolisAnnealing();                                                 // Metropolis algorithm for the parameters space

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

        double m_beta;                                                              // Temperature of the Simulated annealing schema
        double m_energy;                                                            // Energy of the accepted configuration
        double m_currentEnergy, m_currentError;                                     // Energy of the proposed moved (not necessarily accepted) and relative error
        int m_N;                                                                    // Number of blocks
        int m_L;                                                                    // Number of throws
        TrialWaveFunction *m_f;                                                     // Trial wave function object
        Random *m_rnd;                                                              // Pointer to the random class
        DoubleDwellPotential *m_d;                                                  // Potential object
};


#endif // __Annealing__