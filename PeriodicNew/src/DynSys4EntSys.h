#pragma once
#include "EstEntSG.h"

using namespace std;
using namespace arma;

class DynSys4EntSys {
public:
	int dim;
	virtual vec f(const vec &cx) = 0;
	virtual vec CoordTrans(const vec &ox) =0;
	virtual mat Amat(const vec &cx) = 0;	  // Df(cx)
	DynSys4EntSys(int _dim) : dim(_dim) {};
};


// concrete systems


class UnitCircle : public DynSys4EntSys {
public:
    UnitCircle(void) : DynSys4EntSys(3) {};
    vec f(const vec &cx);
    vec CoordTrans(const vec &ox);
    //mat Vmat(const vec &cx);
    mat Amat(const vec &cx);
};

class Orbit3D : public DynSys4EntSys {
public:
    Orbit3D(void) : DynSys4EntSys(3) {};
    vec f(const vec &cx);
    vec CoordTrans(const vec &ox);
    //mat Vmat(const vec &cx);
    mat Amat(const vec &cx);
};
class SprottSystem : public DynSys4EntSys {
public:
    SprottSystem(void) : DynSys4EntSys(3) {};
    vec f(const vec &cx);
    vec CoordTrans(const vec &ox);
    //mat Vmat(const vec &cx);
    mat Amat(const vec &cx);
};
class repressilator : public DynSys4EntSys {
public:
    repressilator(void) : DynSys4EntSys(3) {};
    vec f(const vec &cx);
    vec CoordTrans(const vec &ox);
    //mat Vmat(const vec &cx);
    mat Amat(const vec &cx);
};
