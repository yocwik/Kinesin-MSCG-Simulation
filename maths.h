//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// This file contains mathematical structures like vectors and
// and matrices and their operations.
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
#ifndef MATHS_INCLUDED
#define MATHS_INCLUDED

#include <cmath>
#include <cstdlib>
#include <cstdio>

#define PI 3.14159265

using namespace std;

//////////////////////////////////////////////////////////////
///////////////////MATHEMATICAL STRUCTURES////////////////////

//defining a data structure for a vector in 3D space
struct mvector {
    double x,y,z;
    mvector();                              //initialize vector
    mvector(double,double,double);          //initialize vector with specific values
    double length();                        //calculate vector's length
    double square();                        //calculate vector's squared length
    double& operator() (int);               //overloading the () operator to access vector's members
};

//defining a data structure for a 3 by 3 matrix
struct matrix {
    double e [3][3];
    matrix();                               //initialize matrix
    void reset();                           //reset matrix values to zero
    double& operator()(int,int);            //overloading the () operator to access the matrix members
};

//defining a data structure for a n by m matrix
struct gmatrix {
    double ** e;
    int m,n;
    gmatrix(int,int);                       //initialize gmatrix
    ~gmatrix();                             //deconstruction of gmatrix
    void reset();                           //reset gmatrix values to zero
    double& operator()(int,int);            //overloading the () operator to access the gmatrix members
};

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////VECTOR OPERATIONS///////////////////////

//initializing mvector
mvector::mvector()
{
    x = 0;
    y = 0;
    z = 0;
}

//initializing values
mvector::mvector(double a,double b,double c)
{
    x = a;
    y = b;
    z = c;
}

//calculate vector length
double mvector::length()
{
    return sqrt(x*x+y*y+z*z);
}

//calculate vector square length
double mvector::square()
{
    return (x*x+y*y+z*z);
}

//defining index access to vector components
double& mvector::operator()(int i)
{
    switch (i)
    {
        case 0:
            return x;
            break;
        case 1:
            return y;
            break;
        case 2:
            return z;
            break;
        default:
            abort();
    }
}

//defining vector addition
mvector operator+(mvector u,mvector v)
{
    mvector res;
    res.x = u.x + v.x;
    res.y = u.y + v.y;
    res.z = u.z + v.z;

    return res;
}

void operator+=(mvector &u,mvector v)
{
    u.x += v.x;
    u.y += v.y;
    u.z += v.z;
}

//defining vector multiplication by a scalar
mvector operator*(mvector u,double p)
{
    mvector res;
    res.x = u.x*p;
    res.y = u.y*p;
    res.z = u.z*p;

    return res;
}

mvector operator*(double p,mvector u)
{
    mvector res;
    res.x = u.x*p;
    res.y = u.y*p;
    res.z = u.z*p;

    return res;
}

void operator*=(mvector &u,double p)
{
    u.x *= p;
    u.y *= p;
    u.z *= p;
}

//defining vector subtraction
mvector operator-(mvector u,mvector v)
{
    mvector res;
    res.x = u.x - v.x;
    res.y = u.y - v.y;
    res.z = u.z - v.z;

    return res;
}

void operator-=(mvector &u,mvector v)
{
    u.x -= v.x;
    u.y -= v.y;
    u.z -= v.z;
}

//defining vector division by a scalar
mvector operator/(mvector u,double p)
{
    mvector res;
    res.x = u.x/p;
    res.y = u.y/p;
    res.z = u.z/p;

    return res;
}

void operator/=(mvector &u,double p)
{
    u.x /= p;
    u.y /= p;
    u.z /= p;
}

//defining dot product
double operator*(mvector a,mvector b)
{
    return (a.x*b.x + a.y*b.y + a.z*b.z);
}

//defining cross product
mvector operator^(mvector a,mvector b)
{
    mvector r;
    r.x = a.y*b.z - a.z*b.y;
    r.y = a.z*b.x - a.x*b.z;
    r.z = a.x*b.y - a.y*b.x;
    return r;
}

//defining a function to set the values of a vector
mvector set_vector(double i,double j,double k)
{
    mvector r;

    r.x = i;
    r.y = j;
    r.z = k;

    return r;
}

//defining a function to calculate the distance between the end points of two vectors
double dist(mvector a,mvector b)
{
    mvector r;
    r = a-b;
    return r.length();
}

//defining a function to calculate the square of the distance between the end points of two vectors
double dist_sqr(mvector a,mvector b)
{
    mvector r;
    r = a-b;
    return r.square();
}

//defining a function to calculate the unit vector in the same direction as a given vector
mvector unit_vector(mvector a)
{
    mvector r;
    r = a/a.length();
    return r;
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////MATRIX OPERATIONS///////////////////////

//initialize matrix
matrix::matrix()
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            e[i][j] = 0;
        }
    }
}

//reset matrix values to 0
void matrix::reset()
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            e[i][j] = 0;
        }
    }
}

//defining index access to matrix components
double& matrix::operator()(int n,int m)
{
    return e[n][m];
}

//defining matrix addition
matrix operator+(matrix a,matrix b)
{
    matrix res;

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            res(i,j) = a(i,j) + b(i,j);
        }
    }

    return res;
}

void operator+=(matrix &a,matrix b)
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            a(i,j) += b(i,j);
        }
    }
}

//defining matrix multiplication by a scalar
matrix operator*(matrix a,double p)
{
    matrix res;

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            res(i,j) = a(i,j)*p;
        }
    }

    return res;
}

matrix operator*(double p,matrix a)
{
    matrix res;

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            res(i,j) = a(i,j)*p;
        }
    }

    return res;
}

void operator*=(matrix &a,double p)
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            a(i,j) *= p;
        }
    }
}

//defining matrix subtraction
matrix operator-(matrix a,matrix b)
{
    matrix res;

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            res(i,j) = a(i,j) - b(i,j);
        }
    }

    return res;
}

void operator-=(matrix &a,matrix b)
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            a(i,j) -= b(i,j);
        }
    }
}

//defining matrix division by a scalar
matrix operator/(matrix a,double p)
{
    matrix res;

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            res(i,j) = a(i,j)/p;
        }
    }

    return res;
}

void operator/=(matrix &a,double p)
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            a(i,j) /= p;
        }
    }
}

//defining matrix multiplication
matrix operator*(matrix a,matrix b)
{
    matrix res;

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            res(i,j) = a(i,0)*b(0,j) + a(i,1)*b(1,j) + a(i,2)*b(2,j);
        }
    }

    return res;
}

//defining a function that returns an identity matrix
matrix eye()
{
    matrix I;
    I.reset();

    for(int i=0;i<3;i++)
    {
        I(i,i) = 1;
    }

    return I;
}

//defining a function that return the transpose matrix
matrix transpose(matrix A)
{
    matrix T;

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            T(i,j) = A(j,i);
        }
    }

    return T;
}
/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////GMATRIX OPERATIONS//////////////////////

//initializing gmatrix
gmatrix::gmatrix(int a,int b)
{
    m = a;
    n = b;

    e = new double*[m];
    for (int i=0;i<m;i++)
    {
        e[i] = new double[n];
    }

    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            e[i][j] = 0;
        }
    }
}

//deconstruction of gmatrix
gmatrix::~gmatrix()
{
    for (int i=0;i<m;i++)
    {
        delete [] e[i];
    }

    delete [] e;

    m = 0;
    n = 0;
}

//reset matrix values to 0
void gmatrix::reset()
{
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            e[i][j] = 0;
        }
    }
}

//defining index access to gmatrix components
double& gmatrix::operator()(int a,int b)
{
    return e[a][b];
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
/////////////////MATRIX & VECTOR OPERATIONS///////////////////

//defining vector multiplication by a matrix
mvector operator*(matrix a,mvector b)
{
    mvector res;

    for(int i=0;i<3;i++)
    {
        res(i) = a(i,0)*b(0) + a(i,1)*b(1) + a(i,2)*b(2);
    }

    return res;
}

//defining tensor multiplication
matrix tensor(mvector u,mvector v)
{
    matrix res;

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            res(i,j) = u(i)*v(j);
        }
    }

    return res;
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////ROTATION OPERATIONS/////////////////////

//rotates a vector v around unit vector u by theta radians
mvector rotate_vector_rodrigues(mvector v,mvector u,double theta)
{
    mvector r;

    r = v*cos(theta) + (u^v)*sin(theta) + u*(u*v)*(1-cos(theta));

    return r;
}

//rotates a vector v around a vector u, whose length is the rotation angle theta in radians
mvector rotate_vector_rodrigues(mvector v,mvector u)
{
    mvector r;
    double theta;

    theta = u.length();

    if(theta!=0.0)
    {
        u /= theta;

        r = v*cos(theta) + (u^v)*sin(theta) + u*(u*v)*(1-cos(theta));
    } else {
        r = v;
    }

    return r;
}

//rotation matrix for cartesian rotation vector k
matrix rotation_matrix_rodrigues(mvector k)
{
    matrix A,K;
    double Db1;

    Db1 = k.length();

    if(Db1!=0) k /= Db1;

    K(0,1) = -k.z;
    K(0,2) = k.y;
    K(1,0) = k.z;
    K(1,2) = -k.x;
    K(2,0) = -k.y;
    K(2,1) = k.x;

    A = eye() + sin(Db1)*K + (1-cos(Db1))*K*K;

    return A;
}

//sums consecutive rotation vectors. u and v give rotation axes. th1/2 give angels. The function returns a vector k3 whose
//direction defines the summed rotation axis and whose magnitude is equal to the rotation angel in radians.
//The summation itself is done using Gibbs representation of the rotation vectors.
mvector add_rotation_gibbs (mvector u, mvector v)
{
    mvector k1,k2,k3;
    double Db1,th1,th2;

    th1 = u.length();
    th2 = v.length();

    if(th1!=0.0) u /= th1;
    if(th2!=0.0) v /= th2;

    k1 = tan(th1/2.0)*u;
    k2 = tan(th2/2.0)*v;

    k3 = (k1 + k2 - (k1^k2))/(1-(k1*k2));

    Db1 = 2.0*atan(k3.length());

    k3 = Db1*unit_vector(k3);

    return k3;
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
////////////////////////IO OPERATIONS/////////////////////////

void mvector_print(mvector v)
{
    printf("%8.6lf %8.6lf %8.6lf\n",v(0),v(1),v(2));
}

void matrix_print(matrix m)
{
    printf("%8.6lf %8.6lf %8.6lf\n",m(0,0),m(0,1),m(0,2));
    printf("%8.6lf %8.6lf %8.6lf\n",m(1,0),m(1,1),m(1,2));
    printf("%8.6lf %8.6lf %8.6lf\n",m(2,0),m(2,1),m(2,2));
}

void gmatrix_fprint(gmatrix& gm,string path)
{
    FILE * f1;
    f1 = fopen(path.c_str(),"w");

    for(int i=0;i<gm.m;i++)
    {
        for(int j=0;j<gm.n;j++)
        {
            fprintf(f1,"%3.3lf ",gm(i,j));
        }

        fprintf(f1,"\n");
    }

    fclose(f1);
}

//defining function to read matrix from file
matrix matrix_read_file(string path)
{
    FILE * f1;

    matrix Mt1;

    f1 = fopen(path.c_str(),"r");

    fscanf(f1,"%lf%lf%lf\n",&Mt1(0,0),&Mt1(0,1),&Mt1(0,2));
    fscanf(f1,"%lf%lf%lf\n",&Mt1(1,0),&Mt1(1,1),&Mt1(1,2));
    fscanf(f1,"%lf%lf%lf\n",&Mt1(2,0),&Mt1(2,1),&Mt1(2,2));

    fclose(f1);

    return Mt1;
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

#endif // MATHS_INCLUDED
