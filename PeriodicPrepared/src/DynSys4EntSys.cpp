#include "DynSys4EntSys.h"
#include <armadillo>

using namespace arma;

// Lorenz
Lorenz::Lorenz(double _sigma, double _rho, double _beta) : sigma(_sigma), rho(_rho), beta(_beta), DynSys4EntSys(3) {
	K_Rad = sqrt(beta / 2.0) * (sigma + rho);
}

vec Lorenz::f(const vec &cx) {
	vec fx(dim);
	fx(0) = sigma * (-cx(0) + cx(1));
	fx(1) = rho * cx(0) - cx(1) - cx(0) * cx(2);
	fx(2) = -beta * cx(2) + cx(0) * cx(1);
	return fx;
}

mat Lorenz::Amat(const vec &cx) {
	return { {-sigma, sigma, 0}, {rho - cx(2), -1, -cx(0)}, {cx(1), cx(0), -beta} };
}

vec Lorenz::CoordTrans(const vec &ox) {
	vec sox = ox % vec{ K_Rad, datum::pi, 2.0 * datum::pi };
	vec cx(dim);
	cx(0) = sox(0) * sin(sox(1)) * cos(sox(2));
	cx(1) = sox(0) * sin(sox(1)) * sin(sox(2));
	cx(2) = sox(0) * cos(sox(1)) + (sigma + rho);
	return cx;
}

// BouncingBall
vec BouncingBall::f(const vec &cx) {
	vec fx(dim);
	fx(0) = cx(0) + cx(1);
	fx(1) = gamma * cx(1) - delta * cos(cx(0) + cx(1));
	return fx;
}

mat BouncingBall::Amat(const vec &cx) {
	return { {1.0, 1.0}, {delta * sin(cx(0) + cx(1)), gamma + delta * sin(cx(0) + cx(1))} };
}

vec BouncingBall::CoordTrans(const vec &ox) {
	double x_min = -delta / (1 - gamma), x_max = 2.0 * datum::pi, y_max = delta / (1.0 - gamma);
	vec cx(dim);
	cx(0) = x_min + ox(0) * (x_max - x_min);
	cx(1) = ox(1) * y_max;
	return cx;
}

// Henon a=1.4, b=0.3
Henon::Henon(void) : DynSys4EntSys(2) {
	vec q00 = { -1.862, 1.96 };
	vec q10 = { -1.484, -2.3333 };
	vec q11 = { 1.743, -0.6533 };
	vec q01 = { 1.848, 0.6267 };
	mat Mq = { {q10(0) - q00(0), q01(0) - q00(0)}, {q10(1) - q00(1), q01(1) - q00(1)} };
	vec cd = solve(Mq, q11 - q00);
	TransLambda = cd(0) + cd(1) - 1;
	Transl = vec{ 1 - cd(1), 1 - cd(0) };
	TranssMat = Mq * diagmat(cd) + q00 * Transl.t();
	Transq00 = q00;
	// (TranssMat * x + Translambda * Transq00)/ (dot(Transl, x)  + Translambda)
}

vec Henon::f(const vec &cx) {
	vec fx(dim);
	fx(0) = 1.4 - pow(cx(0), 2) + 0.3 * cx(1);
	fx(1) = cx(0);
	return fx;
}

mat Henon::Amat(const vec &cx) {
	return { {-2.0 * cx(0), 0.3}, {1.0, 0.0} };
}

vec Henon::CoordTrans(const vec &ox) {
	return (TranssMat * ox + TransLambda * Transq00) / (dot(Transl, ox) + TransLambda);
}

//Van der Pol system mu=3
vec VanDerPol::f(const vec &cx) {
    vec fx(dim);
    fx(0) = -cx(1);
    fx(1) = cx(0)-3*cx(1)-3*pow(cx(0),2)*cx(1);
    return fx;
}

mat VanDerPol::Amat(const vec &cx) {
    return { {0, -1}, {1-6*cx(0)*cx(1), -3-3*pow(cx(0),2)} };
}

vec VanDerPol::CoordTrans(const vec &ox) {
    double x_min=-1, y_min=-2.8;
    double a = 2.2, b=-0.2, c =2.8, d=2.8;
    vec cx(dim);
    cx(0) = a*ox(0)+b*ox(1)+x_min;
    cx(1) = c*ox(0)+d*ox(1)+y_min;
    return cx;
}


// UnitCircle
vec UnitCircle::f(const vec &cx) {
    vec fx(dim);
    fx(0) = cx(0)-pow(cx(0),3)-cx(0)*pow(cx(1),2)-cx(1);
    fx(1) = cx(1)-pow(cx(0),2)*cx(1)-pow(cx(1),3)+cx(0);
    return fx;
}

/*
mat UnitCircle::Amat(const vec &cx) {
    return { {1-3*pow(cx(0),2)-pow(cx(1),2), -2*cx(0)*cx(1)-1}, {-2*cx(0)*cx(1)+1, 1-pow(cx(0),2)-3*pow(cx(1),2)} };
}*/

mat UnitCircle::Amat(const vec& cx) {
    vec fx(dim);
    fx(0) = cx(0)-pow(cx(0),3)-cx(0)*pow(cx(1),2)-cx(1);
    fx(1) = cx(1)-pow(cx(0),2)*cx(1)-pow(cx(1),3)+cx(0);
    mat amat;
	
    amat.resize ( 2, 2 ) ;
    amat(0,0) =1-3*pow(cx(0),2)-pow(cx(1),2);
    amat(0,1)= -2*cx(0)*cx(1)-1;
    amat(1,0) = -2*cx(0)*cx(1)+1;
    amat(1,1) = 1-pow(cx(0),2)-3*pow(cx(1),2);
    return amat - (fx*fx.t()*(amat+amat.t()))/pow(norm(fx,2),2);
}

vec UnitCircle::CoordTrans(const vec& ox) {
    double innerRad = 0.59, outerRad = 1.41;
    vec sox = (ox % vec{ outerRad - innerRad, 2 * datum::pi });
    vec cx(dim);
    cx(0) = (innerRad + sox(0)) * cos(sox(1));
    cx(1) = (innerRad + sox(0)) * sin(sox(1));
    return cx;
}

vec Orbit3D::f(const vec& cx) {
	vec fx(dim);
	fx(0) = cx(0) * (1 - pow(cx(0), 2) - pow(cx(1), 2)) - cx(1) + 0.1 * cx(1) * cx(2);
	fx(1) = cx(1) * (1 - pow(cx(0), 2) - pow(cx(1), 2)) + cx(0);
	fx(2) = -cx(2) + cx(0) * cx(1);
	return fx;
}

mat Orbit3D::Amat(const vec& cx) {
	return { {1 - 3 * pow(cx(0), 2) - pow(cx(1), 2), -2 * cx(0) * cx(1) - 1 + 0.1 * cx(2), 0.1 * cx(1)},
			{-2 * cx(0) * cx(1) + 1, 1 - 3 * pow(cx(1), 2) - pow(cx(0), 2), 0},
			{cx(1), cx(0), -1} };
}

vec Orbit3D::CoordTrans(const vec& ox) {
	double r = 0.07, R = 1.0;
	double wonkiness = 0.12;
	vec sox = (ox % vec{ r, 2 * datum::pi, 2 * datum::pi });
	vec cx(dim);
	cx(0) = (R + sox(0) * cos(sox(1))) * cos(sox(2));
	cx(1) = (R + sox(0) * cos(sox(1))) * sin(sox(2));
	cx(2) = sox(0) * sin(sox(1)) - wonkiness * cos(2 * sox(2));
	return cx;
}

vec Repressilator::f(const vec& cx) {
	vec fx(dim);
	fx(0) = mu / (1 + pow(cx(1), s)) - cx(0);
	fx(1) = mu / (1 + pow(cx(2), s)) - cx(1);
	fx(2) = mu / (1 + pow(cx(0), s)) - cx(2);
	return fx;
}

mat Repressilator::Amat(const vec& cx) {

	vec fx = f(cx);

	mat amat;
	amat.resize(3, 3);
	amat(0, 0) = -1;
	amat(0, 1) = (-mu * s * pow(cx(1), s - 1)) / (pow((1 + pow(cx(1), s)), 2));
	amat(0, 2) = 0;
	amat(1, 0) = 0;
	amat(1, 1) = -1;
	amat(1, 2) = (-mu * s * pow(cx(2), s - 1)) / (pow((1 + pow(cx(2), s)), 2));
	amat(2, 0) = (-mu * s * pow(cx(0), s - 1)) / (pow((1 + pow(cx(0), s)), 2));
	amat(2, 1) = 0;
	amat(2, 2) = -1;

	return amat - (fx * fx.t() * (amat + amat.t())) / pow(norm(fx, 2), 2);
}

vec Repressilator::CoordTrans(const vec& ox) {
	double scale = mu;
	double shift = 0;
	vec cx = scale * ox + shift * vec{ 1, 1, 1 };
	return cx;
}