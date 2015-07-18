#ifndef EQUATION_H
#define EQUATION_H

#include "math.h"
const double si=40.;
const double h=0.0001;
double exact(const double & x, const double & y)
{
    //return 1.;
    //return sin(M_PI*x)*sin(M_PI*y)+sin(M_PI*x*2.)*sin(M_PI*y*2.);
    return pow(M_E,-si*(pow(x-1.,2)+pow(y-1.,2)));
}

double sigma(const double & x, const double & y)
{
    //return 1.;
    return 1.+pow(M_E,-si*(pow(x-1.,2)+pow(y-1.,2)));
    //return 1.+exact(x,y);
}

double source(const double & x, const double & y)
{
    //return 1.;
    //return (1.+2.*M_PI*M_PI)*sin(M_PI*x)*sin(M_PI*y)+(1.+8.*M_PI*M_PI)*sin(M_PI*x*2.)*sin(M_PI*y*2.);
   //return (4*si+1.-4*pow(si*(x-1.),2)-4*pow(si*(y-1.),2))*exact(x,y);
    //return exact(x,y)+(4*si-4*pow(si*(x-1.),2)-4*pow(si*(y-1.),2))*exact(x,y)*sigma(x,y)
    //       -(4*pow(si*(x-1.),2)+4*pow(si*(y-1.),2))*exact(x,y)*exact(x,y) ;
    return exact(x,y)-(sigma(x+0.5*h,y)*(exact(x+h,y)-exact(x,y))-sigma(x-0.5*h,y)*(exact(x,y)-exact(x-h,y)))/(h*h)
             -(sigma(x,y+0.5*h)*(exact(x,y+h)-exact(x,y))-sigma(x,y-0.5*h)*(exact(x,y)-exact(x,y-h)))/(h*h);
}

double kappa(const double & x, const double & y)
{
    return 1.;
}

#endif
