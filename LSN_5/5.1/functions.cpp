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
double Hydrogen210::Evaluate(arma::vec v) const {
    double theta = acos(v[2]/norm(v,2));
    return pow(m_a0, -5) / ( 32 * M_PI ) * pow(arma::norm(v, 2), 2) * exp( - arma::norm(v,2) / m_a0 ) * pow( cos(theta) , 2);
}
