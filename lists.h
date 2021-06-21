//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// This file contains functions for creating and manipulating
// neighbor lists.
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
#ifndef LISTS_INCLUDED
#define LISTS_INCLUDED

#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>

#include "maths.h"
#include "ranfun.h"

using namespace std;

//////////////////////////////////////////////////////////////
///////////////////////DATA STRUCTURES////////////////////////

//defining a data structure to contain a pair of integer indexes
struct index_pair
{
    int a,b;

    index_pair()
    {
        a = 0;
        b = 0;
    }
};

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
////////////////////SIMULATION PARAMETERS/////////////////////

//neighbor list additional cutoff distance
double RM = 4.0;

//neighbor list cutoff distance
double NLcutoff = 30.0;

//neighbor list total cutoff distance squared
double NLR2  = (NLcutoff + RM)*(NLcutoff + RM);

//neighbor list update condition
bool NLupdate = true;

//neighbor list sum of maximal displacements
double rmax2 = 0;

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
///////////////////////NEIGHBOR LISTS/////////////////////////

//create a record of the displacement of a molecule
vector<mvector> create_displacement_record(int n)
{
    vector<mvector> dr;

    for(int i=0;i<n;i++)
    {
        dr.push_back(set_vector(0,0,0));
    }

    return dr;
}

//reset the displacement record
void reset_diplacement_record(vector<mvector>& pos0,vector<mvector>& pos)
{
    int n = (int) pos.size();

    for(int i=0;i<n;i++)
    {
        pos0[i] = pos[i];
        pos[i] = set_vector(0,0,0);
    }
}

//create a molecular neighbor list
vector<index_pair> create_neighbor_list(const vector<mvector>& pos,const vector<index_pair>& il)
{
    vector<index_pair> nl;
    int n = (int) il.size();

    for(int i=0;i<n;i++)
    {
        if(dist_sqr(pos[il[i].a],pos[il[i].b]) < NLR2)
        {
            nl.push_back(il[i]);
        }
    }

    return nl;
}

//create an inter-molecular neighbor list
vector<index_pair> create_neighbor_list(const vector<mvector>& pos1,const vector<mvector>& pos2,const vector<index_pair>& il)
{
    vector<index_pair> nl;
    int n = (int) il.size();

    for(int i=0;i<n;i++)
    {
        if(dist_sqr(pos1[il[i].a],pos2[il[i].b]) < NLR2)
        {
            nl.push_back(il[i]);
        }
    }

    return nl;
}

//update the displacement list
void update_displacement_list(vector<double>& displacement,const vector<mvector>& pos0,const vector<mvector>& pos)
{
    mvector dr;

    int n = (int) pos.size();

    displacement.clear();

    for(int i=0;i<n;i++)
    {
        dr = pos[i] - pos0[i];

        displacement.push_back(dr.square());
    }
}

//decide whether to update the neighbor list or not
bool update_neighbor_list(vector<double> displacement)
{
    bool update = false;

    //sort displacements
    std::sort(displacement.begin(),displacement.end());
    std::reverse(displacement.begin(),displacement.end());

    //determine maximal displacements
    rmax2 = sqrt(displacement[0]) + sqrt(displacement[1]);

    //determine update status
    if(rmax2 > RM)
    {
        update = true;
    } else {
        update = false;
    }

    return update;
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

#endif // SIMULATION_INCLUDED
