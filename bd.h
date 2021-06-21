//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// This file defines variables and functions for integrating
// Brownian dynamics of point particles
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
#ifndef BD_INCLUDED
#define BD_INCLUDED

#include <string>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <ctime>

#include "maths.h"
#include "ranfun.h"

using namespace std;

//////////////////////////////////////////////////////////////
//////////////////////////VARIABLES///////////////////////////

//time integration step (ps)
double h = 0.1;

//KT (pN*A)
double KT = 41.14;

//viscosity (pN*ps/(A^2))
double eta = 30.0;

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
/////////////////BROWNIAN DYNAMICS INTEGRATORS////////////////

//initiating random number seed
long init = time(NULL);

//defining a function to perform a single integration step following the standard Brownian dynamics algorithm
void integrate_brownian(vector<mvector>& beads,const vector<mvector>& f,const vector<mvector>& fluc,const vector<double> fc,const vector<double> rc,int first,int last)
{
    for(int i=first;i<last;i++)
    {
        beads[i] += fc[i]*f[i] + rc[i]*fluc[i];
    }
}

//defining a function to perform a single integration step following the standard Brownian dynamics algorithm
void integrate_brownian_BAOAB(vector<mvector>& beads,const vector<mvector>& f,const vector<mvector>& fluc1,const vector<mvector>& fluc2,const vector<double> fc,const vector<double> rc,int first,int last)
{
    for(int i=first;i<last;i++)
    {
        beads[i] += fc[i]*f[i] + rc[i]*0.5*(fluc1[i]+fluc2[i]);
    }
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

#endif
