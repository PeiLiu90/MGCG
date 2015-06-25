#ifndef EQUATION_H
#define EQUATION_H

#include "math.h"

double exact(const double & x, const double & y)
{
    //return x*x+y*y;
    //return sin(M_PI*x)*sin(M_PI*y)+sin(M_PI*x*2.)*sin(M_PI*y*2.);
    return pow(M_E,-20.*(pow(x-1.,2)+pow(y-1.,2)));
}

double source(const double & x, const double & y)
{
    //return 1.;
    //return x*x+y*y;
    //return (1.+2.*M_PI*M_PI)*sin(M_PI*x)*sin(M_PI*y)+(1.+8.*M_PI*M_PI)*sin(M_PI*x*2.)*sin(M_PI*y*2.);
    return (81.-pow(40.*(x-1.),2)-pow(40.*(y-1.),2))*pow(M_E,-20.*(pow(x-1.,2)+pow(y-1.,2)));
}

double kappa(const double & x, const double & y)
{
    //return sqrt(x*x+y*y)*2.;
    return 1.;
}

double sigma(const int & i, const double & x, const double & y)
{
    if(i==0||i==3)
    {
        return 1.;
    }
    else
    {
        return 0.;
    }
}
#endif
