#include "functions.h"
#include "random.h"

using namespace std;

TrialWaveFunction::TrialWaveFunction(Random *rnd){

    /*
    Trial wave function constructor
    */

    m_mu = 1;
    m_sigma = 1.;
    m_rnd = rnd;
    m_position = 0;                                                                  // initial position x=0
    m_stepLength = 2;                                                                // step of the proposed Metropolis move
}

double TrialWaveFunction::Evaluate(double x) const {

    /*
    Evaluation method
    */

    return pow(abs(exp(-pow(x - m_mu, 2) / (2 * pow(m_sigma, 2))) + exp(-pow(x + m_mu, 2) / (2 * pow(m_sigma, 2)))), 2);
}

double TrialWaveFunction::EvaluateNoModulus(double x) const {

    /*
    Evaluation without modulus
    */

    return abs(exp(-pow(x - m_mu, 2) / (2 * pow(m_sigma, 2))) + exp(-pow(x + m_mu, 2) / (2 * pow(m_sigma, 2))));
}

double TrialWaveFunction::SecondDerivative(double x) const {

    /*
    Implementation of the second derivative used in the variational Montecarlo code for the Trial wave function
    */

    double minusExp = exp(-0.5 * (pow(x-m_mu,2)/(pow(m_sigma,2))));
    double plusExp = exp(-0.5 * (pow(x + m_mu, 2) / (pow(m_sigma, 2))));

    return ((-1 / pow(m_sigma, 2)) * minusExp) + ((-1 / pow(m_sigma, 2)) * plusExp) + ((pow(x - m_mu, 2) / pow(m_sigma, 4)) * minusExp) + ((pow(x + m_mu, 2) / pow(m_sigma, 4)) * plusExp);
}

void TrialWaveFunction::MetropolisUniform() {

    /*
    Performs Metropolis update dirrectly with the trial wave function
    */

    double tempPosition = m_position + m_rnd->Rannyu(-1, 1) * m_stepLength;
    double alpha = min(1., (Evaluate(tempPosition) / Evaluate(m_position)));

    // accepting the configuration with probability \alpha:
    double p = m_rnd->Rannyu();
    if (p < alpha){
        m_position = tempPosition;
    }
}

void TrialWaveFunction::Equilibrate(int nblocks, int L) {

    /*
    Performs the equilibration routine needed for Metropolis
    */
    
    for (int j = 0; j < nblocks; j++)
    {
        for (int i = 0; i < L; i++)
        {
            MetropolisUniform();
        }
    }
}

//=======================================================================================================================//

DoubleDwellPotential::DoubleDwellPotential(){

    /*
    Constructor
    */

    m_a = -2.5;
    m_b = 1;
}

double DoubleDwellPotential::Evaluate(double x) const {

    /*
    Evaluation method
    */

    return m_a * pow(x, 2) + m_b * pow(x,4);

}

void DoubleDwellPotential::SetCoeffs(arma::vec coeff){
    m_a = coeff[0];
    m_b = coeff[1];
}
