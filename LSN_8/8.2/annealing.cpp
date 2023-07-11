#include "annealing.h"

using namespace std;

Annealing::Annealing(TrialWaveFunction *TWF, Random *rnd, DoubleDwellPotential *DDP) {

    /*
    Constructor of the Annealing class
    Inputs:
    - TrialWaveFunction *TWF: pointer to a trial wave function object, used to access all its methods.
    - Random *rnd: pointer to the Random class
    - DoubleDwellPotential *DDP: pointer to the function object which represents the potential I am using
    */

    m_beta = 1;
    m_energy = 1;
    m_currentEnergy = m_energy;
    m_currentError = 0;
    m_f = TWF;
    m_rnd = rnd;
    m_d = DDP;
    m_N = 0;
    m_L = 0;
}

void  Annealing::MetropolisAnnealing() {

    /*
    Carries out the Metropolis evolution for the simulated annealing (Not the pdf sampling)
    */

    double oldMu = m_f->Get_Mu();                                                   // Save the old parameters that will be changed if the Metropolis step is admissable
    double oldSigma = m_f->Get_Sigma();

    m_f->Set_Mu(abs(oldMu + m_rnd->Rannyu(-1, 1) * .5 * (1/m_beta)));               // Update the parameters of the trial wave function: this is the proposed move
    m_f->Set_Sigma(abs(oldSigma + m_rnd->Rannyu(-1, 1) * .25 * (1/m_beta)));

    // variables for the blocking average
    double runningSum = 0.;
    double runningSquared = 0.;
    double error = 0.;
    double integral = 0.;

    m_f->Equilibrate(100, pow(10, 3));                                              // Equilibrate the system for the given trial wave function

    for (int i = 0; i < m_N; i++) {                                                 // Loop over the number of blocks

        integral = 0;

        for (int j = 0; j < m_L; j++){                                              // Loop over the number of throws in each block

            m_f->MetropolisUniform();                                               // updates my position according to my trial wave function distribution and the Metropolis evolution
            integral += ((-0.5 * m_f->SecondDerivative(m_f->Get_Position())) / m_f->EvaluateNoModulus(m_f->Get_Position())) + m_d->Evaluate(m_f->Get_Position());

        }

        // Data - blocking, I/O managment
        integral /= m_L;
        runningSum = ((runningSum * i) + integral) / (i + 1);                       // at every iteration I have to re-multiplicate the previous division
        runningSquared = ((runningSquared * i) + pow(integral, 2)) / (i + 1);
        error = Error(runningSum, runningSquared, i);

    }

    m_currentEnergy = runningSum;
    m_currentError = error;

    double alpha = min(1., (exp(- m_beta * (m_currentEnergy - m_energy))));         // Check wheter the proposed moved is acceptable through Metropolis

    // Accepting the configuration with probability \alpha (this means rejecting it with probability (1-\alpha)):
    double p = m_rnd->Rannyu();

    if (alpha > p){                                                                 // Hold the computed configuration

        m_energy = runningSum;
    }

    else                                                                            //Get back to old configuration
    {
        m_f ->Set_Mu(oldMu);
        m_f-> Set_Sigma(oldSigma);
    }
}

