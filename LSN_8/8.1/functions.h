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
        virtual double Evaluate (double x) const =0;                                  // By defining it virtual, we make sure that every child class needs to implement this method
        virtual double Evaluate (arma::vec v) const = 0;
};

class Cosine: public Functions {

    /*
    Derived from Functions, it builds a Cosine obj.
    */

    private:
        double m_a;
        double m_b;
        double m_c;

    public:
        Cosine();                                                                     // empty constructor
        Cosine(double a, double b, double c);                                         // Constructor
        ~Cosine() {;};                                                                //empty destructor

        virtual double Evaluate (double x) const {return m_a * cos(m_b * x + m_c);}   // Implementation of the virtual method

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
        Parabola();                                                                   // empty constructor
        Parabola(double a, double b, double c);                                       // Constructor
        ~Parabola() { ; };                                                            // empty destructor

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
        HydrogenGS();                                                                 // empty constructor
        ~HydrogenGS() { ; };                                                          // empty destructor

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
        Hydrogen210();                                                                // empty constructor
        ~Hydrogen210() { ; };                                                         // empty destructor

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

class TrialWaveFunction: public Functions {

    /*
    Derived from functions, it build a trial wave function obj used for the Variational montecarlo code.
    */

    public: 
        TrialWaveFunction();                                                          // empty constructor
        ~TrialWaveFunction() {;};                                                     // empty destructor

        virtual double Evaluate(arma::vec v) const
        {
            cout << "This WF is for a 1D model, hence it works with doubles not vectors. There is something wrong" << endl;
            return 0;
        };
        virtual double Evaluate(double x) const;                                      // Evaluation method
        double EvaluateNoModulus(double x) const;                                     // Evaluation without modulus
        double SecondDerivative(double x) const;                                      // Second derivative
        void Set_bohr_radius(double a) { m_a0 = a; }
        void Set_Mu(double mu) {m_mu = mu;}
        void Set_Sigma (double sigma) {m_sigma = sigma;}
        double Get_bohr_radius() { return m_a0; }
        double Get_Mu() {return m_mu;}
        double Get_Sigma() {return m_sigma;}

    private:
        double m_a0, m_mu, m_sigma;

};

class DoubleDwellPotential: public Functions {

    /*
    Derived from FUnctions, it build a dwell-shape potential obj.
    */

    public:
        DoubleDwellPotential();                                                       // empty constructor
        ~DoubleDwellPotential() {;};                                                  // empty destructor

        virtual double Evaluate(arma::vec v) const{
            cout << "The potential is 1D, please give a 1D INPUT, not a VECTOR" << endl;
            return 0;
        };
        virtual double Evaluate(double x) const;                                      // Evaluation
        void SetCoeffs(arma::vec coeff);

    private:

        double m_a, m_b;

};

#endif // __Functions__