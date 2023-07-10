#include "functions.h"

using namespace std;

Cosine::Cosine()
{

    /*
    Cosine empty constructor
    */

    m_a = 1;
    m_b = 1;
    m_c = 0;
}

Cosine::Cosine(double a, double b, double c)
{

    /*
    Cosine constructor
    */

    m_a = a;
    m_b = b;
    m_c = c;
}

Parabola::Parabola()
{

    /*
    Parabola empty constructor
    */

    m_a = 1;
    m_b = 1;
    m_c = 1;
}

Parabola::Parabola(double a, double b, double c)
{

    /*
    Parabola constructor
    */

    m_a = a;
    m_b = b;
    m_c = c;
}

HydrogenGS::HydrogenGS()
{

    /*
    Hydrogen GS constructor
    */

    m_a0 = 0.0529 * pow(10, -9);
}

Hydrogen210::Hydrogen210()
{

    /*
    Hydrogen excited state constuctor
    */

    m_a0 = 0.0529 * pow(10, -9);
}

//Trial wave function constructor
TrialWaveFunction::TrialWaveFunction(){
    m_a0 = 0.0529 * pow(10, -9);
    m_mu = 0.;
    m_sigma = 1.;
}

double Hydrogen210::Evaluate(arma::vec v) const {

    /*
    Evaluation method for the excited state of the Hydrogen pdf.
    */

    double theta = acos(v[2]/norm(v,2));
    return pow(m_a0, -5) / ( 32 * M_PI ) * pow(arma::norm(v, 2), 2) * exp( - arma::norm(v,2) / m_a0 ) * pow( cos(theta) , 2);
}

double TrialWaveFunction::Evaluate(double x) const {

    /*
    Evaluation method for the trial wafe function pdf.
    */

    return pow(abs(exp(-pow(x - m_mu, 2) / (2 * pow(m_sigma, 2))) + exp(-pow(x + m_mu, 2) / (2 * pow(m_sigma, 2)))), 2);
}

double TrialWaveFunction::EvaluateNoModulus(double x) const {

    /*
    This is the evalutaion without the modulus
    */

    return abs(exp(-pow(x - m_mu, 2) / (2 * pow(m_sigma, 2))) + exp(-pow(x + m_mu, 2) / (2 * pow(m_sigma, 2))));
}

double TrialWaveFunction::SecondDerivative(double x) const{

    /*
    Implementation of the second derivative used in the variational Montecarlo code for the Trial wave function
    */

    double minusExp = exp(-0.5 * (pow(x - m_mu, 2) / (pow(m_sigma, 2))));
    double plusExp = exp(-0.5 * (pow(x + m_mu, 2) / (pow(m_sigma, 2))));

    return ((-1 / pow(m_sigma, 2)) * minusExp) + ((-1 / pow(m_sigma, 2)) * plusExp) + ((pow(x - m_mu, 2) / pow(m_sigma, 4)) * minusExp) + ((pow(x + m_mu, 2) / pow(m_sigma, 4)) * plusExp);
}

DoubleDwellPotential::DoubleDwellPotential(){

    /*
    Constructor of the dwell-shaped potential
    */

    m_a = -2.5;
    m_b = 1;
}

double DoubleDwellPotential::Evaluate(double x) const {

    /*
    Evaluation method for the dwell-shaped potential
    */

    return m_a * pow(x, 2) + m_b * pow(x,4);

}

void DoubleDwellPotential::SetCoeffs(arma::vec coeff){

    /*
    Simple set method
    */

    m_a = coeff[0];
    m_b = coeff[1];
}
