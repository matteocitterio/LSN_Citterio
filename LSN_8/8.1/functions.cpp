#include "functions.h"

using namespace std;

//Cosine constructors
Cosine::Cosine(){
    m_a = 1;
    m_b = 1;
    m_c = 0;
}

Cosine::Cosine(double a, double b, double c){
    m_a = a;
    m_b = b;
    m_c = c;
}

// Parabola constructors
Parabola::Parabola()
{
    m_a = 1;
    m_b = 1;
    m_c = 1;
}

Parabola::Parabola(double a, double b, double c)
{
    m_a = a;
    m_b = b;
    m_c = c;
}

//Hydrogen Ground State constructors
HydrogenGS::HydrogenGS()
{
    m_a0 = 0.0529 * pow(10, -9);
}
Hydrogen210::Hydrogen210()
{
    m_a0 = 0.0529 * pow(10, -9);
}

//Trial wave function constructor
TrialWaveFunction::TrialWaveFunction(){
    m_a0 = 0.0529 * pow(10, -9);
    m_mu = 0.;
    m_sigma = 1.;
}

//Evaluate function
double Hydrogen210::Evaluate(arma::vec v) const {
    double theta = acos(v[2]/norm(v,2));
    return pow(m_a0, -5) / ( 32 * M_PI ) * pow(arma::norm(v, 2), 2) * exp( - arma::norm(v,2) / m_a0 ) * pow( cos(theta) , 2);
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

    double minusExp = exp(-0.5 * (pow(x - m_mu, 2) / (pow(m_sigma, 2))));
    double plusExp = exp(-0.5 * (pow(x + m_mu, 2) / (pow(m_sigma, 2))));

    return ((-1 / pow(m_sigma, 2)) * minusExp) + ((-1 / pow(m_sigma, 2)) * plusExp) + ((pow(x - m_mu, 2) / pow(m_sigma, 4)) * minusExp) + ((pow(x + m_mu, 2) / pow(m_sigma, 4)) * plusExp);
}

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
