//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// This file defines functions for manipulating coarse grained
// representation of rigid objects.
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
#ifndef RBDCG_INCLUDED
#define RBDCG_INCLUDED

#include <string>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>

#include "maths.h"
#include "ranfun.h"

using namespace std;

//////////////////////////////////////////////////////////////
///////////////////////////DATA I/O///////////////////////////

//defining a function to read the positions of cg beads from file
void read_cg_beads(vector<mvector>& beads,string path)
{
    mvector Mv1;

    FILE * f1;

    f1 = fopen(path.c_str(),"r");

    beads.clear();

    while(!feof(f1))
    {
        fscanf(f1,"%lf %lf %lf\n",&Mv1.x,&Mv1.y,&Mv1.z);
        beads.push_back(Mv1);
    }

    fclose(f1);
}

//defining a function to read the radii of cg beads from file
void read_cg_radii(vector<double>& rad,string path)
{
    double Db1;

    FILE * f1;

    f1 = fopen(path.c_str(),"r");

    rad.clear();

    while(!feof(f1))
    {
        fscanf(f1,"%lf\n",&Db1);
        rad.push_back(Db1);
    }

    fclose(f1);
}

//defining a function to write the positions of cg beads in XYZ format to file
void write_cg_beads_xyz(vector<mvector>& beads,string path)
{
    int In1;

    In1 = (int) beads.size();

    FILE * f1;

    f1 = fopen(path.c_str(),"w");

    fprintf(f1,"%d\n\n",In1);

    for(int i=0;i<In1;i++)
    {
        fprintf(f1,"C %lf %lf %lf\n",beads[i].x,beads[i].y,beads[i].z);
    }

    fclose(f1);
}

//defining a function to write a frame of cg bead positions in XYZ format to a file
void write_cg_beads_xyz_frame(vector<mvector>& beads,string path)
{
    int In1;

    In1 = (int) beads.size();

    FILE * f1;

    f1 = fopen(path.c_str(),"a");

    fprintf(f1,"%d\n\n",In1);

    for(int i=0;i<In1;i++)
    {
        fprintf(f1,"C %lf %lf %lf\n",beads[i].x,beads[i].y,beads[i].z);
    }

    fclose(f1);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
////////////////////MECHANICAL OPERATIONS/////////////////////

//defining a function to update the bead positions of a rigid object
void update_rbd_coordinates(vector<mvector>& coord,vector<mvector>& coord_init,mvector pos,mvector rot)
{
    int In1;

    In1 = (int) coord.size();

    for(int i=0;i<In1;i++)
    {
        coord[i] = rotate_vector_rodrigues(coord_init[i],rot) + pos;
    }
}

//defining a function to update the bead positions of a rigid object within an index range
void update_rbd_coordinates(vector<mvector>& coord,vector<mvector>& coord_init,mvector pos,mvector rot,int a, int b)
{
    for(int i=a;i<b;i++)
    {
        coord[i] = rotate_vector_rodrigues(coord_init[i],rot) + pos;
    }
}

//defining a function to calculate the resultant force
mvector calc_RB_resultant_force(const vector<mvector>& f,int a,int b)
{
    mvector Mv1;

    for(int i=a;i<b;i++)
    {
        Mv1 += f[i];
    }

    return Mv1;
}

//defining a function to calculate the resultant torque around the origin
mvector calc_RB_resultant_torque(const vector<mvector>& t,const vector<mvector>& r,const mvector& ori,int a,int b)
{
    mvector Mv1;

    for(int i=a;i<b;i++)
    {
        Mv1 += (r[i]-ori)^t[i];
    }

    return Mv1;
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
/////////////////////////////MISC/////////////////////////////

//defining a function to create a vmd topology representation file for a cg molecule
void create_VMD_topo_CG(string path,vector<double> radii,int a)
{
    int In1;

    In1 = (int) radii.size();


    FILE * f1;

    f1 = fopen(path.c_str(),"w");

    fprintf(f1,"mol top %d\n\n",a);

    fprintf(f1,"mol modstyle 0 %d VDW 1.0\n\n",a);

    for(int i=0;i<In1;i++)
    {
        fprintf(f1,"set sel [atomselect top \"index %d\"]\n",i);
        fprintf(f1,"$sel set radius {%lf}\n",radii[i]);
    }

    fclose(f1);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

#endif
