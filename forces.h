//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// This file contains variables and functions for calculating
// and applying forces.
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
#ifndef FORCES_INCLUDED
#define FORCES_INCLUDED

#include <iostream>
#include <cmath>

#include "maths.h"

using namespace std;

//////////////////////////////////////////////////////////////
//////////////////////////CONSTANTS///////////////////////////

//defining boltzman's constant (KCal/(MOL*K))
const double KB = 1.98720431e-3;

//defining avogadro's number
const double MOL = 6.02214179e+23;

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
/////////////////////////CONVERSIONS//////////////////////////

//Note: All conversions assume that distance is in angstroms

//Energy
//number of J in one KCal/Mol
const double KCALPMOL_TO_JOULE = 4184/MOL;

//number of KCal/Mol in one J
const double JOULE_TO_KCALPMOL = MOL/4184;

//Force
//number of pN in one KCal/(Mol*A)
const double KCALPMOLA_TO_PN = 4184e+22/MOL;

//number of KCal/(Mol*A) in one pN
const double PN_TO_KCALPMOLA = MOL/4184e+22;

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////FORCE PARAMETERS////////////////////////

//fene stiffness
double k_fene;

//fene typical distance
double R0_fene;

//fene force coefficients
double fene1,fene2;

//repulsion energy constant
double eps_r;

//repulsion force coefficient
double rep1;

//////////////////////////////////////////////////////////////
///////////////////////INITIALIZATION/////////////////////////

//defining a function to obtain the force parameters
void set_force_param()
{
    k_fene = 20.0*KCALPMOLA_TO_PN; //convert Kcal per mol to pN units

    R0_fene = 2.0; //A

    eps_r = 1.0*KCALPMOLA_TO_PN; //convert Kcal per mol to pN units
}

//defining a function to calculate the force coefficients
void set_force_coeff()
{
    //constant calculations
    fene1 = k_fene*R0_fene*R0_fene;
    fene2 = R0_fene*R0_fene;
    rep1 = 6*eps_r;
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
/////////////////////////FORCE DATA///////////////////////////

//defining a function to reset forces
void reset_forces(vector<mvector>& forces)
{
    for(int i=0;i<(int) forces.size();i++)
    {
        forces[i] = set_vector(0,0,0);
    }
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
///////////////////////////FORCES/////////////////////////////

//calculate harmonic force
mvector calc_harmonic(mvector bead1,mvector bead2,double k,double x0)
{
    mvector r;
    double dr,x,f;

    //calculate unit vector connecting atoms a and b
    r = bead1 - bead2;
    dr = r.length();
    r /= dr;

    //calculating force coordinate
    x = dr - x0;

    //calculate force magnitude
    f = -k*x;

    //finalize force
    r *= f;

    return r;
}

//calculate fene force
mvector calc_fene(mvector bead1,mvector bead2,double x0)
{
    mvector r;
    double dr,x,f;

    //calculate unit vector connecting atoms a and b
    r = bead1 - bead2;
    dr = r.length();
    r /= dr;

    //calculating force coordinate
    x = dr - x0;

    //calculate force magnitude
    f = -fene1*x/(fene2-x*x);

    //finalize force
    r *= f;

    return r;
}

//calculate modified repulsive force
mvector calc_repulsion(mvector bead1,mvector bead2,double rad1,double rad2,double sig)
{
    mvector r;
    double dr,f,shft;

    //determine force length scale
    if(rad1+rad2<=sig)
    {
         //calculate unit vector connecting atoms a and b
        r = bead1 - bead2;
        dr = r.square();

        //calculate force magnitude
        f = rep1*pow(sig,6)/pow(dr,4);

        //finalize force
        r *= f;
    } else {
        //calculate unit vector connecting atoms a and b
        r = bead1 - bead2;
        dr = r.length();
        r /= dr;

        //calculating shift
        shft = rad1 + rad2 - sig;

        //calculate force magnitude
        f = rep1*pow(sig,6)/pow(dr-shft,7);

        //finalize force
        r *= f;
    }

    return r;
}

//calculate modified repulsive force
mvector calc_repulsion_MT(mvector bead1,double rad1,double rad_MT,double sig)
{
    mvector r;
    double dr,f,shft;

    //calculate unit vector connecting atoms a and b
    r = bead1;
    r.x = 0;
    dr = r.length();
    r /= dr;

    //calculating shift
    shft = rad1 + rad_MT - 2.0*sig;

    //calculate force magnitude
    f = rep1*pow(sig,6)/pow(dr-shft,7);

    //finalize force
    r *= f;

    return r;
}

//calculate angular potential force
void calc_angular(mvector& MvA,mvector& MvB,mvector& MvfA,mvector& MvfB,double ang0,double k)
{
    double DbA,DbB,DbC,ang,coeff;
    mvector MvC;

    DbA = MvA.length();
    DbB = MvB.length();

    ang = acos(MvA*MvB/(DbA*DbB));

    coeff = -k*(ang-ang0);

    MvC = MvB^MvA;

    DbC = MvC.length();

    if(DbC==0.0)
    {
        MvfA = set_vector(0,0,0);
        MvfB = set_vector(0,0,0);
    } else {
        MvfA = coeff*(MvC^MvA)/(DbA*DbA*DbC);
        MvfB = coeff*(MvB^MvC)/(DbB*DbB*DbC);
    }
}

//calculate dihedral potential force
void calc_dihedral(mvector& MvA,mvector& MvB,mvector& MvC,mvector& Mvf1,mvector& Mvf2,mvector& Mvf3,mvector& Mvf4,double ang0,double k,double& ang_rec,double& shft)
{
    double Db1,Db2,Db3,ang,coeff,Dbt1,Dbt2;
    mvector Mv1,Mv2,MvfA,MvfB,MvfC;

    Mv1 = MvA^MvC;
    Mv2 = MvB^MvC;

    Db1 = Mv1.length();
    Db2 = Mv2.length();
    Db3 = MvC.length();

    Dbt1 = (Mv1*Mv2)/(Db1*Db2);

    Dbt2 = ((Mv2^Mv1)*MvC)/(Db1*Db2*Db3);

    if(abs(Dbt1)>1.0) Dbt1 /= abs(Dbt1);

    if(abs(Dbt2)>1.0) Dbt2 /= abs(Dbt2);

    ang = acos(Dbt1);

    if(asin(Dbt2)<0.0) ang *= -1.0;

    if(ang+shft-ang_rec>PI)
    {
        shft -= 2.0*PI;
    }

    if(ang+shft-ang_rec<-PI)
    {
        shft += 2.0*PI;
    }

    ang += shft;

    ang_rec = ang;

    //E = 0.5*k*(ang-ang0)^2
    coeff = -k*(ang-ang0);

    MvfA = -(Db3/(Db1*Db1))*Mv1;
    MvfB = (Db3/(Db2*Db2))*Mv2;
    MvfC = ((MvA*MvC)/(Db1*Db1*Db3))*Mv1 - ((MvB*MvC)/(Db2*Db2*Db3))*Mv2;

    Mvf1 = coeff*MvfA;
    Mvf2 = coeff*(MvfC - MvfA);
    Mvf3 = -coeff*(MvfC + MvfB);
    Mvf4 = coeff*MvfB;
}

//calculate attraction well force
mvector calc_attraction(mvector bead1,mvector bead2,double k,double u)
{
    mvector r;
    double dr,f;

    //calculate unit vector connecting atoms a and b
    r = bead1 - bead2;
    dr = r.square();

    //calculate force magnitude
    f = -2.0*k*u/pow(1.0+k*dr,2.0);

    //finalize force
    r *= f;

    return r;
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

#endif // FORCES_INCLUDED
