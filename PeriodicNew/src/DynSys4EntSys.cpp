#include "DynSys4EntSys.h"



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

mat UnitCircle::Amat(const vec &cx){

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



// 3D example periodic orbit
vec Orbit3D::f(const vec &cx) {
    vec fx(dim);
    fx(0) = cx(0)*(1-pow(cx(0),2)-pow(cx(1),2))-cx(1)-0.1*cx(1)*cx(2);
    fx(1) = cx(1)*(1-pow(cx(0),2)-pow(cx(1),2))+cx(0);
    fx(2) = -cx(2)+cx(0)*cx(1);
    return fx;
}


mat Orbit3D::Amat(const vec &cx){

    vec fx(dim);
    fx(0) = cx(0)*(1-pow(cx(0),2)-pow(cx(1),2))-cx(1)-0.1*cx(1)*cx(2);
    fx(1) = cx(1)*(1-pow(cx(0),2)-pow(cx(1),2))+cx(0);
    fx(2) = -cx(2)+cx(0)*cx(1);

    mat amat;

    amat.resize ( 3, 3 ) ;
    amat(0,0) = 1-3*pow(cx(0),2)-pow(cx(1),2);
    amat(0,1) = -2*cx(0)*cx(1)-1+0.1*cx(2);
    amat(0,2) = 0.1*cx(1);
    amat(1,0) = 1-pow(cx(0),2)-3*pow(cx(1),2);
    amat(1,1) = -2*cx(0)*cx(1);
    amat(1,2) = 0;
    amat(2,0) = cx(1);
    amat(2,1) = cx(0);
    amat(2,2) = -1;
    return amat - (fx*fx.t()*(amat+amat.t()))/pow(norm(fx,2),2);
}

/*
 // Area from Iman, works well
vec Orbit3D::CoordTrans(const vec& ox) {
    double r = 0.20, R = 1.0;
    double wonkiness = 0.20;
    vec sox = (ox % vec{ r, 2 * datum::pi, 2 * datum::pi });
    vec cx(dim);
    cx(0) = (R + sox(0) * cos(sox(1))) * cos(sox(2));
    cx(1) = (R + sox(0) * cos(sox(1))) * sin(sox(2));
    cx(2) = sox(0) * sin(sox(1)) - wonkiness * cos(2 * sox(2));
    return cx;
}*/

//Based on calculations from Peter
// Need to test these parameters for more iterarions
vec Orbit3D::CoordTrans(const vec& ox) {
    double r = 0.66, R = sqrt(1.33);
    double z = 1.84;
    vec sox = (ox % vec{ r, 2 * datum::pi, 2 * datum::pi });
    vec cx(dim);
    cx(0) = (R + sox(0) * cos(sox(1))) * cos(sox(2));
    cx(1) = (R + sox(0) * cos(sox(1))) * sin(sox(2));
    cx(2) = sox(0) * sin(sox(1))*z/r;
    return cx;
}

// Sprott system
// Not yet successful
vec SprottSystem::f(const vec &cx) {
    double mu = 1;
    vec fx(dim);
    fx(0) = pow(cx(1),2)-cx(2)-mu*cx(0);
    fx(1) = pow(cx(2),2)-cx(0)-mu*cx(1);
    fx(2) = pow(cx(0),2)-cx(1)-mu*cx(2);
    return fx;
}


mat SprottSystem::Amat(const vec &cx){

    double mu = 1;
    vec fx(dim);
    fx(0) = pow(cx(1),2)-cx(2)-mu*cx(0);
    fx(1) = pow(cx(2),2)-cx(0)-mu*cx(1);
    fx(2) = pow(cx(0),2)-cx(1)-mu*cx(2);

    mat amat;

    amat.resize ( 3, 3 ) ;
    amat(0,0) = -mu;
    amat(0,1) = 2*cx(1);
    amat(0,2) = -1;
    amat(1,0) = -1;
    amat(1,1) = -mu;
    amat(1,2) = 2*cx(2);
    amat(2,0) = 2*cx(0);
    amat(2,1) = -1;
    amat(2,2) = -mu;
    return amat - (fx*fx.t()*(amat+amat.t()))/pow(norm(fx,2),2);
}

vec SprottSystem::CoordTrans(const vec& ox) {
    double mu =1;
    vec sox = ox % vec{ 0.8, datum::pi, 2.0 * datum::pi };
    vec cx(dim);
    cx(0) = sox(0) * sin(sox(1)) * cos(sox(2))*0.5;
    cx(1) = sox(0) * sin(sox(1)) * sin(sox(2))*0.5;
    cx(2) = sox(0) * cos(sox(1));
    return cx;
}



// repressilator system
// See other project: EntEstSG_periodic_repressilator
vec repressilator::f(const vec &cx) {
    double mu = 1.0;
    double s  = 20.0;
    vec fx(dim);
    fx(0) = mu/(1+pow(cx(1),s))-cx(0);
    fx(1) = mu/(1+pow(cx(2),s))-cx(1);
    fx(2) = mu/(1+pow(cx(0),s))-cx(2);
    return fx;
}


mat repressilator::Amat(const vec &cx){

    double mu = 1.0;
    double s  = 20.0;
    vec fx(dim);
    fx(0) = mu/(1+pow(cx(1),s))-cx(0);
    fx(1) = mu/(1+pow(cx(2),s))-cx(1);
    fx(2) = mu/(1+pow(cx(0),s))-cx(2);

    mat amat;

    amat.resize ( 3, 3 ) ;
    amat(0,0) = -1;
    amat(0,1) = (-mu*s*pow(cx(1),s-1))/(pow((1+pow(cx(1),s)),2));
    amat(0,2) = 0;
    amat(1,0) = 0;
    amat(1,1) = -1;
    amat(1,2) = (-mu*s*pow(cx(2),s-1))/(pow((1+pow(cx(2),s)),2));
    amat(2,0) = (-mu*s*pow(cx(0),s-1))/(pow((1+pow(cx(0),s)),2));
    amat(2,1) = 0;
    amat(2,2) = -1;
    return amat - (fx*fx.t()*(amat+amat.t()))/pow(norm(fx,2),2);
}

/*vec repressilator::CoordTrans(const vec& ox) {
    double mu =1;
    vec sox = ox % vec{ 0.8, datum::pi, 2.0 * datum::pi };
    vec cx(dim);
    cx(0) = sox(0) * sin(sox(1)) * cos(sox(2))*0.5;
    cx(1) = sox(0) * sin(sox(1)) * sin(sox(2))*0.5;
    cx(2) = sox(0) * cos(sox(1));
    return cx;
}*/
vec repressilator::CoordTrans(const vec& ox) {
    double mu = 1.0;
    double s = 20.0;
    double scale = mu - mu / (2 + pow(mu, s));
    double shift = mu / (2 + pow(mu, s));
    vec cx = scale * ox + shift * vec{ 1, 1, 1 };
    return cx;
}

