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

HydrogenGS::HydrogenGS() {

    /*
    Hydrogen GS constructor
    */

    m_a0 = 0.0529 * pow(10, -9);
}

Hydrogen210::Hydrogen210() {

    /*
    Hydrogen excited state constuctor
    */

    m_a0 = 0.0529 * pow(10, -9);
}

double Hydrogen210::Evaluate(arma::vec v) const {

    /*
    Evaluation method for the excited state of the Hydrogen pdf.
    */

    double theta = acos(v[2]/norm(v,2));
    return pow(m_a0, -5) / ( 32 * M_PI ) * pow(arma::norm(v, 2), 2) * exp( - arma::norm(v,2) / m_a0 ) * pow( cos(theta) , 2);
}
