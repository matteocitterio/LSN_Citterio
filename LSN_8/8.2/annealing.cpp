#include "annealing.h"

using namespace std;

Annealing::Annealing(TrialWaveFunction *TWF, Random *rnd, DoubleDwellPotential *DDP){
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

void  Annealing::MetropolisAnnealing(){

    // Save the old parameters
    double oldMu = m_f->Get_Mu();
    double oldSigma = m_f->Get_Sigma();

    // Update the parameters:
    m_f->Set_Mu(abs(oldMu + m_rnd->Rannyu(-1, 1) * .5 * (1/m_beta)));
    m_f->Set_Sigma(abs(oldSigma + m_rnd->Rannyu(-1, 1) * .25 * (1/m_beta)));

    //cout << "Trying new config: MU=" << m_f->Get_Mu() << " SIGMA= " << m_f->Get_Sigma() << endl;

    // variables for the blocking average
    double runningSum = 0.;
    double runningSquared = 0.;
    double error = 0.;
    double integral = 0.;

    m_f->Equilibrate(100, pow(10, 3));

    for (int i = 0; i < m_N; i++){
        integral = 0;

        for (int j = 0; j < m_L; j++){

            // updates my position according to my trial wave function distribution
            m_f->MetropolisUniform();
            // Now we need to compute the energy
            integral += ((-0.5 * m_f->SecondDerivative(m_f->Get_Position())) / m_f->EvaluateNoModulus(m_f->Get_Position())) + m_d->Evaluate(m_f->Get_Position());
        }
        integral /= m_L;
        runningSum = ((runningSum * i) + integral) / (i + 1); // at every iteration I have to re-multiplicate the previous division
        runningSquared = ((runningSquared * i) + pow(integral, 2)) / (i + 1);
        error = Error(runningSum, runningSquared, i);
        //cout << "Running sum:" << runningSum<< endl;
    }

    m_currentEnergy = runningSum;
    m_currentError = error;

    //cout << "OLD ENERGY=" << m_energy << " NEW ENERGY=" << m_currentEnergy << endl;
    //cout << " exponential " << exp(-1 * m_beta * (m_currentEnergy - m_energy)) << endl;

    //Now we perform Metropolis:
    double alpha = min(1., (exp(- m_beta * (m_currentEnergy - m_energy))));

    //cout << "ALPHA: " << alpha << endl;
    // accepting the configuration with probability \alpha (this means rejecting it with probability (1-\alpha)):
    double p = m_rnd->Rannyu();

    if (alpha > p){

        m_energy = runningSum;
    }

    else                      //Get back to old configuration
    {
        m_f ->Set_Mu(oldMu);
        m_f-> Set_Sigma(oldSigma);
    }
}

