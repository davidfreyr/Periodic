//#pragma once
#include "EstEntSG.h"
#include <armadillo>

using namespace std;
using namespace arma;

class DynSys4EntSys {
public:
	int dim;
    bool isRepressilator = false;
	virtual vec f(const vec &cx) = 0;
	virtual vec CoordTrans(const vec &ox) = 0;
	virtual mat Amat(const vec &cx) = 0;	  // Df(cx)
    double eq = 0;
    vec eqPnt;
	DynSys4EntSys(int _dim) : dim(_dim) {};
};


// concrete systems

class Lorenz : public DynSys4EntSys {
public:
	double sigma, rho, beta, K_Rad;
	Lorenz(double _sigma = 10.0, double _rho = 28.0, double _beta = 8.0 / 3.0);
	vec f(const vec &cx);
	vec CoordTrans(const vec &ox);
	mat Amat(const vec &cx);
};

class BouncingBall : public DynSys4EntSys {
public:
	double delta, gamma;
	BouncingBall(double _delta = 2.0, double _gamma = 0.1)
		: delta(_delta), gamma(_gamma), DynSys4EntSys(2) {};
	vec f(const vec &cx);
	vec CoordTrans(const vec &ox); 
	mat Amat(const vec &cx);
};

class Henon : public DynSys4EntSys {
public:
	double TransLambda;
	vec Transq00, Transl;
	mat TranssMat;
	Henon(void);
	vec f(const vec &cx);
	vec CoordTrans(const vec &ox);
	mat Amat(const vec &cx);
};

class VanDerPol : public DynSys4EntSys {
public:
    VanDerPol(void)
            : DynSys4EntSys(2) {};
    vec f(const vec &cx);
    vec CoordTrans(const vec &ox);
    mat Amat(const vec &cx);
};

class SpeedControl : public DynSys4EntSys {
public:
    SpeedControl(void)
            : DynSys4EntSys(2) {};
    vec f(const vec &cx);
    vec CoordTrans(const vec &ox);
    mat Amat(const vec &cx);
};

class JetEngine : public DynSys4EntSys {
public:
    JetEngine(void)
            : DynSys4EntSys(2) {};
    vec f(const vec &cx);
    vec CoordTrans(const vec &ox);
    mat Amat(const vec &cx);
};

class three_d_example : public DynSys4EntSys {
public:
    three_d_example(void) : DynSys4EntSys(3) {};
    vec f(const vec &cx);
    vec CoordTrans(const vec &ox);
    mat Amat(const vec &cx);
};

class UnitCircle : public DynSys4EntSys {
public:
    UnitCircle(void) : DynSys4EntSys(2) {};
    vec f(const vec &cx);
    vec CoordTrans(const vec &ox);
    //mat Vmat(const vec &cx);
    mat Amat(const vec &cx);
};

class Orbit3D : public DynSys4EntSys {
public:
    Orbit3D(void) : DynSys4EntSys(3) {};
    vec f(const vec& cx);
    vec CoordTrans(const vec& ox);
    mat Amat(const vec& cx);
};

class Repressilator : public DynSys4EntSys {
public:
    double mu = 5.0;
    double s = 3.0;
    vec f(const vec& cx);
    vec CoordTrans(const vec& ox);
    mat Amat(const vec& cx);
    Repressilator(void) : DynSys4EntSys(3) {
        eq = 0.04;
        eqPnt = vec{ eq, eq, eq };
        for (int i = 0; i < 10; i++) {
            vec feqPnt = f(eqPnt);
            mat J;
            J.resize(3, 3);
            J(0, 0) = -1;
            J(0, 1) = (-mu * s * pow(eqPnt(1), s - 1)) / (pow((1 + pow(eqPnt(1), s)), 2));
            J(0, 2) = 0;
            J(1, 0) = 0;
            J(1, 1) = -1;
            J(1, 2) = (-mu * s * pow(eqPnt(2), s - 1)) / (pow((1 + pow(eqPnt(2), s)), 2));
            J(2, 0) = (-mu * s * pow(eqPnt(0), s - 1)) / (pow((1 + pow(eqPnt(0), s)), 2));
            J(2, 1) = 0;
            J(2, 2) = -1;
            eqPnt -= inv(J) * feqPnt;
            cout << eqPnt.t();
        }
        cout << "Equilibrium at " << eqPnt.t();
        eq = eqPnt(0);
    }
};