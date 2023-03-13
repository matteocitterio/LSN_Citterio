#ifndef __Functions__
#define __Functions__

#include <cmath>

using namespace std;

//Parent class used for every basic function
class Functions {

    public:
        virtual double Evaluate (double x) const =0;        //by defining it virtual, we make sure that every child class needs to implement this method

};

//class cosine
class Cosine: public Functions {
    private:
        double m_a;
        double m_b;
        double m_c;

    public:
        Cosine();      //empty constructor
        Cosine(double a, double b, double c);
        ~Cosine() {;};  //empty destructor

        virtual double Evaluate (double x) const {return m_a * cos(m_b * x + m_c);}    //implementation of the virtual method

        void SetA(double a) { m_a = a; }
        void SetB(double b) { m_b = b; }
        void SetC(double c) { m_c = c; }
        double GetA() { return m_a; }
        double GetB() { return m_b; }
        double GetC() { return m_c; }
};

// class parabola
class Parabola : public Functions
{
    private:
        double m_a;
        double m_b;
        double m_c;

    public:
        Parabola(); // empty constructor
        Parabola(double a, double b, double c);
        ~Parabola() { ; }; // empty destructor

        virtual double Evaluate(double x) const { return m_a * pow(x,2) + m_b * x + m_c; } // implementation of the virtual method

        void SetA(double a) { m_a = a; }
        void SetB(double b) { m_b = b; }
        void SetC(double c) { m_c = c; }
        double GetA() { return m_a; }
        double GetB() { return m_b; }
        double GetC() { return m_c; }
};

#endif // __Functions__