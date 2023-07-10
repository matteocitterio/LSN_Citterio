#ifndef __Functions__
#define __Functions__

#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;


class Functions {

    /*
    Parent class used for every basic function
    */

    public:
        virtual double Evaluate (double x) const =0;                                  // by defining it virtual, we make sure that every child class needs to implement this method
        virtual double Evaluate (arma::vec v) const = 0;
};

//class cosine
class Cosine: public Functions {

    /*
    Derived from Functions, it builds a Cosine obj.
    */

    private:
        double m_a;
        double m_b;
        double m_c;

    public:
        Cosine();                                                                      // Empty constructor
        Cosine(double a, double b, double c);                                          // Constructor
        ~Cosine() {;};                                                                 // Empty destructor

        virtual double Evaluate (double x) const {return m_a * cos(m_b * x + m_c);}    // Implementation of the virtual method

        void SetA(double a) { m_a = a; }
        void SetB(double b) { m_b = b; }
        void SetC(double c) { m_c = c; }
        double GetA() { return m_a; }
        double GetB() { return m_b; }
        double GetC() { return m_c; }
};

class Parabola : public Functions {

    /*
    Derived from Functions, it builds a Parabola obj.
    */

    private:
        double m_a;
        double m_b;
        double m_c;

    public:
        Parabola();                                                                     // Empty constructor
        Parabola(double a, double b, double c);                                         // Constructor
        ~Parabola() { ; };                                                              // Empty destructor

        virtual double Evaluate(double x) const { return m_a * pow(x,2) + m_b * x + m_c; } // implementation of the virtual method

        void SetA(double a) { m_a = a; }
        void SetB(double b) { m_b = b; }
        void SetC(double c) { m_c = c; }
        double GetA() { return m_a; }
        double GetB() { return m_b; }
        double GetC() { return m_c; }
};

class HydrogenGS: public Functions {

    /*
    Derived from Functions, it builds a Hydrogen Ground State obj.
    */

    public:
        HydrogenGS();                                                                   // empty constructor
        ~HydrogenGS() { ; };                                                            // empty destructor

        virtual double Evaluate(arma::vec v) const { return pow(m_a0, -3) * exp(-2 * (arma::norm(v,2)) / m_a0) / M_PI; }
        virtual double Evaluate(double x) const
        {
            cout << "This function needs a vector, not a double. There is something wrong" << endl;
            return 0;
        };

        void Set_bohr_radius(double a) { m_a0 = a; }
        double Get_bohr_radius() { return m_a0; }

    private:
        double m_a0;
};

class Hydrogen210: public Functions {

    /*
    Derived from Functions, it builds a Hydrogen Excited State obj.
    */

    public:
        Hydrogen210();                                                                  // empty constructor
        ~Hydrogen210() { ; };                                                           // empty destructor

        virtual double Evaluate(arma::vec v) const;
        virtual double Evaluate(double v) const {
            cout << "This function needs a vector, not a double. There is something wrong" << endl;
            return 0;
        };
        void Set_bohr_radius(double a) { m_a0 = a; }
        double Get_bohr_radius() { return m_a0; }

    private:
        double m_a0;
};

#endif // __Functions__