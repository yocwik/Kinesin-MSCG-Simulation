//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// This file contains variables and functions for building
// the different auxiliary components of the kinesin motor
// system.
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
#ifndef KINSTR_INCLUDED
#define KINSTR_INCLUDED

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

using namespace std;

//////////////////////////////////////////////////////////////
//////////////////////////VARIABLES///////////////////////////

//////////////////////////Coiled Coil/////////////////////////

//number of coiled-coil cg beads
int CC_N_beads;

//cg beads of coiled coil
vector<mvector> CC_beads;

//cg bead radii of coiled coil
double CC_radii;

//coiled coil bead forces
vector<mvector> CC_force;

//coiled coil bead random fluctuations
vector<mvector> CC_fluc1;
vector<mvector> CC_fluc2;

//coiled coil bead force integration coefficients
vector<double> CC_int_fc;

//coiled coil bead random integration coefficients
vector<double> CC_int_rc;

//CC-NL-Junction
//CC-NL Junction cg beads reference positions
vector<mvector> CCNLJ_beads_init;

//CC-NL Junction cg beads
vector<mvector> CCNLJ_beads;

//CC-NL Junction bead forces
vector<mvector> CCNLJ_force;

//CC-NL Junction rigid body frame
RBobj CCNLJ_RB;

//CC-NL Junction rigid body forces
mvector CCNLJ_RB_force;

//CC-NL Junction rigid body torques
mvector CCNLJ_RB_torque;

//CC-NL Junction rigid body random forces
mvector CCNLJ_RB_rand_f1;
mvector CCNLJ_RB_rand_f2;

//CC-NL Junction rigid body random torques
mvector CCNLJ_RB_rand_t1;
mvector CCNLJ_RB_rand_t2;

//CC-NL Junction diffusion tensor
double CCNLJ_DT[6][6];

//CC-NL Junction diffusion tensor Cholesky decomposition
double CCNLJ_DT_Chol[6][6];

//////////////////////////Microtubule/////////////////////////

//MT TBS bead reference positions
vector<mvector> MT_TBS_pos_beads_ref;

//MT TBS rigid object reference position
mvector MT_TBS_pos_ref;

//MT radius to protofilament
double MT_rad = 110.588;

//MT radius of protofilament
double MT_proto_rad = 20.0;

//angular shift between protofilaments
double MT_theta = 2.0*PI/13.0;

//distance between two consecutive alpha tubulin units along the same protofilament
double MT_aa_gap = 81.283;

//distance between the adjacent alpha and beta tubulin units along the same protofilament
double MT_ab_gap = 40.739;

//position shift along the MT axis associated with two adjacent protofilaments
double MT_proto_shift = -9.607;

//number of tubulin repeats
int MT_Nrep = 11;

//////////////////////////Cargo///////////////////////////////

//cargo cg beads reference positions
vector<mvector> CAR_beads_init;

//cargo cg beads
vector<mvector> CAR_beads;

//cargo rigid body
RBobj CAR_RB;

//cargo rigid body forces
mvector CAR_RB_force;

//cargo rigid body torques
mvector CAR_RB_torque;

//kinesin rigid body random forces
mvector CAR_RB_rand_f1;
mvector CAR_RB_rand_f2;

//kinesin rigid body random torques
mvector CAR_RB_rand_t1;
mvector CAR_RB_rand_t2;

//cargo diffusion tensor
double CAR_DT[6][6];

//cargo diffusion tensor Cholesky decomposition
double CAR_DT_Chol[6][6];

//cargo bead forces
vector<mvector> CAR_force;

//cargo radius
double CAR_rad;

//observed angle tracking
double CAR_obs_ang = 0.0;
double CAR_obs_ang_shift = 0.0;

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
///////////////////////INITIALIZATION/////////////////////////

//defining a function to initialize the coiled coil parameters and objects
void CC_initialize_objects()
{
    //number of cg beads (1 is subtracted to account for the CCNLJ)
    CC_N_beads = 30-1;

    //CC bead radius
    CC_radii = 5;

    //initialize forces and fluctuations
    CC_force.assign(CC_N_beads,set_vector(0,0,0));
    CC_fluc1.assign(CC_N_beads,set_vector(0,0,0));
    CC_fluc2.assign(CC_N_beads,set_vector(0,0,0));

    //initialize integration coefficients
    CC_int_fc.assign(CC_N_beads,0.0);
    CC_int_rc.assign(CC_N_beads,0.0);

    for(int i=0;i<CC_N_beads;i++)
    {
        CC_int_fc[i] = h/(6*PI*eta*CC_radii);
        CC_int_rc[i] = sqrt(2.0*h*KT/(6*PI*eta*CC_radii));
    }

    //CC-NL-Junction
    //creating CC-NL-Junction beads (main CC bead, the two ends of the NLs, and the virtual position of the next CC bead)
    CCNLJ_beads_init.assign(4,set_vector(0,0,0));
    CCNLJ_beads.assign(4,set_vector(0,0,0));

    //initializing CC-NL-Junction bead forces
    CCNLJ_force.assign(4,set_vector(0,0,0));

    //calculating diffusion tensor and Cholesky decomposition
    for(int i=0;i<6;i++)
    {
        for(int j=0;j<6;j++)
        {
            if(i==j)
            {
                if(i<3)
                {
                    CCNLJ_DT[i][j] = KT/(6.0*PI*eta*CC_radii);  //A^2/ps
                } else {
                    CCNLJ_DT[i][j] = KT/(8.0*PI*eta*pow(CC_radii,3.0)); //ps^-1
                }

                CCNLJ_DT_Chol[i][j] = sqrt(CCNLJ_DT[i][j]);
            } else {
                CCNLJ_DT[i][j] = 0.0;
                CCNLJ_DT_Chol[i][j] = 0.0;
            }
        }
    }
}

//defining a function to create a coiled coil structure
void CC_create_CC_structure()
{
    mvector MvOri,Mv1;

    //position of the first CC bead (Junction)
    MvOri = set_vector(111.279+40.64,14.8916,147.997+20.0);

    //create CC beads and assign them positions
    CC_beads.assign(CC_N_beads,set_vector(0,0,0));

    for(int i=0;i<CC_N_beads;i++)
    {
        CC_beads[i] = MvOri + (i+1)*set_vector(0,0,2.0*CC_radii);
    }

    //CC-NL-Junction
    //CC-NL-Junction initial position and orientation
    CCNLJ_RB.pos = MvOri;
    CCNLJ_RB.rot = set_vector(0,0,0);

    //position of main CC bead
    CCNLJ_beads_init[0] = set_vector(0,0,0);

    //positions of the NL ends
    CCNLJ_beads_init[1] = set_vector(-2.0,0,-6.6);
    CCNLJ_beads_init[2] = set_vector(2.0,0,-6.6);

    //position of the next CC bead (virtual)
    CCNLJ_beads_init[3] = set_vector(0,0,2.0*CC_radii);

    //update cargo bead coordinates
    update_rbd_coordinates(CCNLJ_beads,CCNLJ_beads_init,CCNLJ_RB.pos,CCNLJ_RB.rot);
}

//defining a function to initialize the MT parameters
//note: has to be used after the KIN parameters have been initialized first but before the TH position is modified from a file
void MT_initialize_objects(vector<mvector> KIN_ref)
{
    //initializing the bead TBS reference positions
    MT_TBS_pos_beads_ref.assign(6,set_vector(0,0,0));

    //setting the bead TBS reference positions
    MT_TBS_pos_beads_ref[0] = KIN_ref[18];
    MT_TBS_pos_beads_ref[1] = KIN_ref[19];
    MT_TBS_pos_beads_ref[2] = KIN_ref[20];
    MT_TBS_pos_beads_ref[3] = KIN_ref[22];
    MT_TBS_pos_beads_ref[4] = KIN_ref[23];
    MT_TBS_pos_beads_ref[5] = KIN_ref[24];

    //setting the rigid object TBS reference position
    MT_TBS_pos_ref = set_vector(111.279,14.8916,147.997);
}

//defining a function to determine the corresponding alpha-beta tubulin dimer for a particular motor position
//note: the function returns an index pair with the first index indicating rotation around the MT and the second translation along the MT
index_pair MT_get_TBS_index(mvector mot)
{
    int In1,In2;
    double ang,pos;
    index_pair MT_ind;

    ang = atan2(mot.y,mot.z) + (MT_theta/2.0);

    In1 = int(floor(ang/MT_theta));

    pos = mot.x - MT_TBS_pos_ref.x - double(In1)*MT_proto_shift + (MT_aa_gap/2.0);

    In2 = int(floor(pos/MT_aa_gap));

    MT_ind.a = In1;
    MT_ind.b = In2;

    return MT_ind;
}

//defining a function to update the TBS location based on an index pair
void MT_update_TBS_pos(vector<mvector>& TBS_beads,index_pair MT_ind)
{
    for(int i=0;i<6;i++)
    {
        TBS_beads[i] = rotate_vector_rodrigues(MT_TBS_pos_beads_ref[i],set_vector(1,0,0),-double(MT_ind.a)*MT_theta) + double(MT_ind.a)*set_vector(MT_proto_shift,0,0) + double(MT_ind.b)*set_vector(MT_aa_gap,0,0);
    }
}

//defining a function that returns the position of the TBS in terms of the rigid body center of motion
mvector MT_get_TBS_RB_pos(index_pair MT_ind)
{
    mvector Mv1;

    Mv1 = rotate_vector_rodrigues(MT_TBS_pos_ref,set_vector(1,0,0),-double(MT_ind.a)*MT_theta) + double(MT_ind.a)*set_vector(MT_proto_shift,0,0) + double(MT_ind.b)*set_vector(MT_aa_gap,0,0);

    return Mv1;
}

//defining a function to create a MT structure file
void MT_create_xyz_file()
{
    //tubulin beads
    vector<mvector> MT_str;

    //create MT
    MT_str.assign(2*MT_Nrep*13,set_vector(0,0,0));

    for(int i=0;i<MT_Nrep;i++)
    {
        for(int j=0;j<13;j++)
        {
            MT_str[i*13+j] = MT_rad*rotate_vector_rodrigues(set_vector(0,0,1),set_vector(1,0,0),MT_theta*(j-6)) - (j-6)*MT_proto_shift*set_vector(1,0,0) + (i-2)*MT_aa_gap*set_vector(1,0,0);
            MT_str[i*13+j+MT_Nrep*13] = MT_rad*rotate_vector_rodrigues(set_vector(0,0,1),set_vector(1,0,0),MT_theta*(j-6)) - (j-6)*MT_proto_shift*set_vector(1,0,0) + (i-2)*MT_aa_gap*set_vector(1,0,0) + MT_ab_gap*set_vector(1,0,0);
        }
    }

    //saving MT structure
    write_cg_beads_xyz(MT_str,"KIN-MT.xyz");
}

//defining a function to create a MT structure file with modified radius
void MT_create_xyz_file_mod_rad(double rmod)
{
    //tubulin beads
    vector<mvector> MT_str;

    //create MT
    MT_str.assign(2*MT_Nrep*13,set_vector(0,0,0));

    for(int i=0;i<MT_Nrep;i++)
    {
        for(int j=0;j<13;j++)
        {
            MT_str[i*13+j] = (MT_rad+rmod)*rotate_vector_rodrigues(set_vector(0,0,1),set_vector(1,0,0),MT_theta*(j-6)) - (j-6)*MT_proto_shift*set_vector(1,0,0) + (i-2)*MT_aa_gap*set_vector(1,0,0);
            MT_str[i*13+j+MT_Nrep*13] = (MT_rad+rmod)*rotate_vector_rodrigues(set_vector(0,0,1),set_vector(1,0,0),MT_theta*(j-6)) - (j-6)*MT_proto_shift*set_vector(1,0,0) + (i-2)*MT_aa_gap*set_vector(1,0,0) + MT_ab_gap*set_vector(1,0,0);
        }
    }

    //saving MT structure
    write_cg_beads_xyz(MT_str,"KIN-MT.xyz");
}

//defining a function to initialize the cargo parameters and objects
void CAR_initialize_objects()
{
    //cargo radius (A)
    CAR_rad = 5000.0;

    //creating cargo beads (main cargo and connection point with CC)
    CAR_beads_init.assign(2,set_vector(0,0,0));
    CAR_beads.assign(2,set_vector(0,0,0));

    //initializing cargo bead forces
    CAR_force.assign(2,set_vector(0,0,0));

    //calculating diffusion tensor and Cholesky decomposition
    for(int i=0;i<6;i++)
    {
        for(int j=0;j<6;j++)
        {
            if(i==j)
            {
                if(i<3)
                {
                    CAR_DT[i][j] = KT/(6.0*PI*eta*CAR_rad);  //A^2/ps
                } else {
                    CAR_DT[i][j] = KT/(8.0*PI*eta*pow(CAR_rad,3.0)); //ps^-1
                }

                CAR_DT_Chol[i][j] = sqrt(CAR_DT[i][j]);
            } else {
                CAR_DT[i][j] = 0.0;
                CAR_DT_Chol[i][j] = 0.0;
            }
        }
    }
}

//defining a function to create the cargo
//note: has to be used after CC has been created
void CAR_create_CAR_structure()
{
    CAR_RB.pos = CC_beads.back() + set_vector(0,0,2.0*CC_radii + CAR_rad);
    CAR_RB.rot = set_vector(0,0,0);

    //positions of main cargo bead
    CAR_beads_init[0] = set_vector(0,0,0);

    //positions of CC connector bead
    CAR_beads_init[1] = set_vector(0,0,-CAR_rad);

    //update cargo bead coordinates
    update_rbd_coordinates(CAR_beads,CAR_beads_init,CAR_RB.pos,CAR_RB.rot);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

#endif
