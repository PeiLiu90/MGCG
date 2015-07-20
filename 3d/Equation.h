#ifndef EQUATION_H
#define EQUATION_H

#include "math.h"
const double si=40.;
const double h=0.0001;
double exact(const double & x, const double & y, const double & z)
{
    //return 1.;
    //return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
    return pow(M_E,-(si*((x-1.)*(x-1.)+(y-1.)*(y-1.)+(z-1.)*(z-1.))));
}

double sigma(const double & x, const double & y, const double & z)
{
    //return 1.;
    return 1.+pow(M_E,-(si*((x-1.)*(x-1.)+(y-1.)*(y-1.)+(z-1.)*(z-1.))));
}

double kappa(const double & x, const double & y, const double & z)
{
    return 1.;
}

double source(const double & x, const double & y, const double & z)
{
    return kappa(x,y,z)*exact(x,y,z)-(sigma(x+0.5*h,y,z)*(exact(x+h,y,z)-exact(x,y,z))-sigma(x-0.5*h,y,z)*(exact(x,y,z)-exact(x-h,y,z)))/(h*h)
                                    -(sigma(x,y+0.5*h,z)*(exact(x,y+h,z)-exact(x,y,z))-sigma(x,y-0.5*h,z)*(exact(x,y,z)-exact(x,y-h,z)))/(h*h)
                                    -(sigma(x,y,z+0.5*h)*(exact(x,y,z+h)-exact(x,y,z))-sigma(x,y,z-0.5*h)*(exact(x,y,z)-exact(x,y,z-h)))/(h*h);
    //return 1.;
}



#endif
