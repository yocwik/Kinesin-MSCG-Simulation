//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// This file defines variables and functions for the
// calculation and saving of simulation observables.
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
#ifndef OBS_INCLUDED
#define OBS_INCLUDED

#include <string>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>

#include "maths.h"
#include "ranfun.h"
#include "bd.h"
#include "rbd.h"
#include "rbdcg.h"
#include "forces.h"
#include "lists.h"
#include "kinstr.h"
#include "kinesin.h"

using namespace std;

//////////////////////////////////////////////////////////////
//////////////////////////VARIABLES///////////////////////////

//path for saving the data
string path_obs = "obs-data";

//observables initial values

//M1 and M2 position
mvector obs_init_M1_pos;
mvector obs_init_M2_pos;

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
/////////////////////////CALCULATIONS/////////////////////////

//defining a function to write the observables to file.
//This function deletes the existing file.
void obs_write_file()
{
    FILE * f1;

    f1 = fopen(path_obs.c_str(),"w");

    //M1 pos
    fprintf(f1,"%lf %lf %lf ",KIN_M1_RB.pos.x,KIN_M1_RB.pos.y,KIN_M1_RB.pos.z);

    //M1 rot
    fprintf(f1,"%lf %lf %lf ",KIN_M1_RB.rot.x,KIN_M1_RB.rot.y,KIN_M1_RB.rot.z);

    //M2 pos
    fprintf(f1,"%lf %lf %lf ",KIN_M2_RB.pos.x,KIN_M2_RB.pos.y,KIN_M2_RB.pos.z);

    //M2 rot
    fprintf(f1,"%lf %lf %lf ",KIN_M2_RB.rot.x,KIN_M2_RB.rot.y,KIN_M2_RB.rot.z);

    //CC-NL-Junction pos
    fprintf(f1,"%lf %lf %lf ",CCNLJ_RB.pos.x,CCNLJ_RB.pos.y,CCNLJ_RB.pos.z);

    //CC-NL-Junction rot
    fprintf(f1,"%lf %lf %lf ",CCNLJ_RB.rot.x,CCNLJ_RB.rot.y,CCNLJ_RB.rot.z);

    //CAR pos
    fprintf(f1,"%lf %lf %lf ",CAR_RB.pos.x,CAR_RB.pos.y,CAR_RB.pos.z);

    //CAR rot
    fprintf(f1,"%lf %lf %lf ",CAR_RB.rot.x,CAR_RB.rot.y,CAR_RB.rot.z);

    //CC 5th bead position
    fprintf(f1,"%lf %lf %lf ",CC_beads[4].x,CC_beads[4].y,CC_beads[4].z);

    //CC end point position
    fprintf(f1,"%lf %lf %lf ",CC_beads[CC_N_beads-1].x,CC_beads[CC_N_beads-1].y,CC_beads[CC_N_beads-1].z);

    //M1 NL end point position
    fprintf(f1,"%lf %lf %lf ",KIN_M1_beads[KIN_N_beads-1].x,KIN_M1_beads[KIN_N_beads-1].y,KIN_M1_beads[KIN_N_beads-1].z);

    //M2 NL end point position
    fprintf(f1,"%lf %lf %lf ",KIN_M2_beads[KIN_N_beads-1].x,KIN_M2_beads[KIN_N_beads-1].y,KIN_M2_beads[KIN_N_beads-1].z);

    //NL docking RMSD
    double Db1,Db2;

    Db1 = 0.0;
    Db2 = 0.0;

    for(int i=0;i<12;i++)
    {
        Db1 += (KIN_M1_beads[i+KIN_N_beads-12] - KIN_M1_dock_pos_nl[i]).square();
        Db2 += (KIN_M2_beads[i+KIN_N_beads-12] - KIN_M2_dock_pos_nl[i]).square();
    }

    Db1 = sqrt(Db1/12.0);
    Db2 = sqrt(Db2/12.0);

    //M1 NL docking
    fprintf(f1,"%lf ",Db1);

    //M2 NL docking
    fprintf(f1,"%lf ",Db2);

    //observed angle
    fprintf(f1,"%lf ",CAR_obs_ang);

    //observed angle shift
    fprintf(f1,"%lf\n",CAR_obs_ang_shift);

    fclose(f1);
}

//defining a function to write the observables to file.
//This function appends data to the existing file.
void obs_write_frame()
{
    FILE * f1;

    f1 = fopen(path_obs.c_str(),"a");

    //M1 pos
    fprintf(f1,"%lf %lf %lf ",KIN_M1_RB.pos.x,KIN_M1_RB.pos.y,KIN_M1_RB.pos.z);

    //M1 rot
    fprintf(f1,"%lf %lf %lf ",KIN_M1_RB.rot.x,KIN_M1_RB.rot.y,KIN_M1_RB.rot.z);

    //M2 pos
    fprintf(f1,"%lf %lf %lf ",KIN_M2_RB.pos.x,KIN_M2_RB.pos.y,KIN_M2_RB.pos.z);

    //M2 rot
    fprintf(f1,"%lf %lf %lf ",KIN_M2_RB.rot.x,KIN_M2_RB.rot.y,KIN_M2_RB.rot.z);

    //CC-NL-Junction pos
    fprintf(f1,"%lf %lf %lf ",CCNLJ_RB.pos.x,CCNLJ_RB.pos.y,CCNLJ_RB.pos.z);

    //CC-NL-Junction rot
    fprintf(f1,"%lf %lf %lf ",CCNLJ_RB.rot.x,CCNLJ_RB.rot.y,CCNLJ_RB.rot.z);

    //CAR pos
    fprintf(f1,"%lf %lf %lf ",CAR_RB.pos.x,CAR_RB.pos.y,CAR_RB.pos.z);

    //CAR rot
    fprintf(f1,"%lf %lf %lf ",CAR_RB.rot.x,CAR_RB.rot.y,CAR_RB.rot.z);

    //CC 5th bead position
    fprintf(f1,"%lf %lf %lf ",CC_beads[4].x,CC_beads[4].y,CC_beads[4].z);

    //CC end point position
    fprintf(f1,"%lf %lf %lf ",CC_beads[CC_N_beads-1].x,CC_beads[CC_N_beads-1].y,CC_beads[CC_N_beads-1].z);

    //M1 NL end point position
    fprintf(f1,"%lf %lf %lf ",KIN_M1_beads[KIN_N_beads-1].x,KIN_M1_beads[KIN_N_beads-1].y,KIN_M1_beads[KIN_N_beads-1].z);

    //M2 NL end point position
    fprintf(f1,"%lf %lf %lf ",KIN_M2_beads[KIN_N_beads-1].x,KIN_M2_beads[KIN_N_beads-1].y,KIN_M2_beads[KIN_N_beads-1].z);

    //NL docking RMSD
    double Db1,Db2;

    Db1 = 0.0;
    Db2 = 0.0;

    for(int i=0;i<12;i++)
    {
        Db1 += (KIN_M1_beads[i+KIN_N_beads-12] - KIN_M1_dock_pos_nl[i]).square();
        Db2 += (KIN_M2_beads[i+KIN_N_beads-12] - KIN_M2_dock_pos_nl[i]).square();
    }

    Db1 = sqrt(Db1/12.0);
    Db2 = sqrt(Db2/12.0);

    //M1 NL docking
    fprintf(f1,"%lf ",Db1);

    //M2 NL docking
    fprintf(f1,"%lf ",Db2);

    //observed angle
    fprintf(f1,"%lf ",CAR_obs_ang);

    //observed angle shift
    fprintf(f1,"%lf\n",CAR_obs_ang_shift);

    fclose(f1);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

#endif
