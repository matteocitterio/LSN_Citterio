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