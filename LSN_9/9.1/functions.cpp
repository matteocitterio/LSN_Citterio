#include "functions.h"
#include "random.h"

using namespace std;

//Trial wave function constructor
TrialWaveFunction::TrialWaveFunction(Random *rnd){

    m_mu = 1;
    m_sigma = 1.;
    m_rnd = rnd;
    m_position = 0;     //initial position x=0
    m_stepLength = 2;   //step of the proposed Metropolis move
}

double TrialWaveFunction::Evaluate(double x) const {
    return pow(abs(exp(-pow(x - m_mu, 2) / (2 * pow(m_sigma, 2))) + exp(-pow(x + m_mu, 2) / (2 * pow(m_sigma, 2)))), 2);
}

//This is the evalutaion without the modulus
double TrialWaveFunction::EvaluateNoModulus(double x) const
{
    return abs(exp(-pow(x - m_mu, 2) / (2 * pow(m_sigma, 2))) + exp(-pow(x + m_mu, 2) / (2 * pow(m_sigma, 2))));
}

double TrialWaveFunction::SecondDerivative(double x) const{

    double minusExp = exp(-0.5 * (pow(x-m_mu,2)/(pow(m_sigma,2))));
    double plusExp = exp(-0.5 * (pow(x + m_mu, 2) / (pow(m_sigma, 2))));

    return ((-1 / pow(m_sigma, 2)) * minusExp) + ((-1 / pow(m_sigma, 2)) * plusExp) + ((pow(x - m_mu, 2) / pow(m_sigma, 4)) * minusExp) + ((pow(x + m_mu, 2) / pow(m_sigma, 4)) * plusExp);
}

void TrialWaveFunction::MetropolisUniform(){

    double tempPosition = m_position + m_rnd->Rannyu(-1, 1) * m_stepLength;
    double alpha = min(1., (Evaluate(tempPosition) / Evaluate(m_position)));

    // accepting the configuration with probability \alpha:
    double p = m_rnd->Rannyu();
    if (p < alpha){
        m_position = tempPosition;
    }
}

void TrialWaveFunction::Equilibrate(int nblocks, int L){
    
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
    m_a = -2.5;
    m_b = 1;
}

double DoubleDwellPotential::Evaluate(double x) const {

    return m_a * pow(x, 2) + m_b * pow(x,4);

}

void DoubleDwellPotential::SetCoeffs(arma::vec coeff){
    m_a = coeff[0];
    m_b = coeff[1];
}