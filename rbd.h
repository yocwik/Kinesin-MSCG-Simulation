//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// This file defines variables and functions for integrating
// Brownian dynamics of rigid bodies
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
#ifndef RBD_INCLUDED
#define RBD_INCLUDED

#include <string>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>

#include "maths.h"
#include "ranfun.h"
#include "bd.h"

using namespace std;

//////////////////////////////////////////////////////////////
//////////////////////////VARIABLES///////////////////////////

//integration coefficient for force/torque
double RBDint_coeff1;

//integration coefficient for random component
double RBDint_coeff2;

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
///////////////////////INITIALIZATION/////////////////////////

//defining a function to initialize the rigid body dynamics parameters
void RBD_initialize_param()
{
    RBDint_coeff1 = h/KT;
    RBDint_coeff2 = sqrt(2.0*h);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
///////////////////////DATA STRUCTURES////////////////////////

//defining a rigid body object data structure
struct RBobj {
    //internal variables

    //position vector
    mvector pos;
    //rotation vector
    mvector rot;

    //operations

    //initialization
    RBobj();
    //member access
    double& operator() (int);               //overloading the () operator to access vector's members
};

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
/////////////////DATA STRUCTURES OPERATIONS///////////////////

//initializing a RBobj
RBobj::RBobj()
{
    pos = set_vector(0,0,0);
    rot = set_vector(0,0,0);
}

//defining index access to rigid body object components
double& RBobj::operator()(int i)
{
    switch (i)
    {
        case 0:
            return pos.x;
            break;
        case 1:
            return pos.y;
            break;
        case 2:
            return pos.z;
            break;
        case 3:
            return rot.x;
            break;
        case 4:
            return rot.y;
            break;
        case 5:
            return rot.z;
            break;
        default:
            abort();
    }
}
/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
////////////////////DIFFUSION TENSOR DATA/////////////////////

//defining a function to read a 6x6 diffusion tensor/decomposition from file and store it in an array
//units: A^2/ps (diffusion tensor)
void read_diff_ten (string path,double (&dt)[6][6])
{
    FILE * f1;

    f1 = fopen(path.c_str(),"r");

    for(int i=0;i<6;i++)
    {
        for(int j=0;j<6;j++)
        {
            fscanf(f1,"%lf",&dt[i][j]);
        }

        fscanf(f1,"\n");
    }

    fclose(f1);
}

//defining a function to multiply the diffusion tensor/decomposition by a scalar
void modify_diff_ten(double (&dt)[6][6],double c)
{
    for(int i=0;i<6;i++)
    {
        for(int j=0;j<6;j++)
        {
            dt[i][j] *= c;
        }
    }
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
////////////////RIGID BODY DYNAMICS OPERATIONS////////////////

//defining a function to generate an mvector of normal random numbers
void get_rand_mvector(mvector& rv,long &init)
{
    for(int i=0;i<3;i++)
    {
        rv(i) = nrand(init);
    }
}

//defining a function to generate an mvector of normal random numbers
mvector get_rand_mvector(long &init)
{
    mvector Mv1;

    for(int i=0;i<3;i++)
    {
        Mv1(i) = nrand(init);
    }

    return Mv1;
}

//defining a function to multiply a 6d vector by the 6x6 tensor from the left. To be used for coupling the different forces
void correlate_force_torque(mvector& f,mvector& t,mvector ref_rot,double (&dt)[6][6])
{
    double ftv0[6],ftv[6];

    //transform to object coordinate system
    f = rotate_vector_rodrigues(f,-1.0*ref_rot);
    t = rotate_vector_rodrigues(t,-1.0*ref_rot);

    for(int i=0;i<3;i++)
    {
        ftv0[i] = f(i);
        ftv0[i+3] = t(i);
    }

    for(int i=0;i<6;i++)
    {
        ftv[i] = 0.0;

        for(int j=0;j<6;j++)
        {
            ftv[i] += dt[i][j]*ftv0[j];
        }
    }

    for(int i=0;i<3;i++)
    {
        f(i) = ftv[i];
        t(i) = ftv[i+3];
    }

    //transform back to universe coordinate system
    f = rotate_vector_rodrigues(f,ref_rot);
    t = rotate_vector_rodrigues(t,ref_rot);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
////////////////RIGID BODY DYNAMICS INTEGRATORS///////////////

//defining a function to perform a single integration step following the standard Brownian dynamics algorithm
void RBD_integrate_brownian(RBobj& obj1,const mvector f,const mvector t,const mvector rf,const mvector rt)
{
    obj1.pos += RBDint_coeff1*f + RBDint_coeff2*rf;
    obj1.rot = add_rotation_gibbs(obj1.rot,RBDint_coeff1*t + RBDint_coeff2*rt);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

#endif
