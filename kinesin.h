//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// This file contains variables and functions for representing
// and manipulating the different objects that are specific to
// the kinesin motor system.
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
#ifndef KIN_INCLUDED
#define KIN_INCLUDED

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

using namespace std;

//////////////////////////////////////////////////////////////
//////////////////////////VARIABLES///////////////////////////

//output structure file paths
string path_str_out_M1 = "KIN-M1.xyz";
string path_str_out_M2 = "KIN-M2.xyz";
string path_str_out_CC = "KIN-CC.xyz";
string path_str_out_CCNLJ = "KIN-CCNLJ.xyz";
string path_str_out_CAR = "KIN-CAR.xyz";

//snapshot input paths
string init_conf_path = "init-conf.dat";

//snapshot output paths
string snapshot_path = "snapshot-conf.dat";

//output VMD topology file path
string path_topo = "topo";

//number of kinesin cg beads
int KIN_N_beads;

//cg beads of kinesin initial position
vector<mvector> KIN_beads_init;

//cg beads of kinesin
vector<mvector> KIN_M1_beads;
vector<mvector> KIN_M2_beads;

//cg docking binding pseudo-beads initial positions
vector<mvector> KIN_dock_pos_cb_init;
vector<mvector> KIN_dock_pos_nl_init;

//cg docking binding pseudo-beads
vector<mvector> KIN_M1_dock_pos_cb;
vector<mvector> KIN_M1_dock_pos_nl;
vector<mvector> KIN_M2_dock_pos_cb;
vector<mvector> KIN_M2_dock_pos_nl;

//cg TBS binding pseudo beads
vector<mvector> MT_TBS_M1_pos;
vector<mvector> MT_TBS_M2_pos;

//TBS indexes
index_pair MT_TBS_M1_ind;
index_pair MT_TBS_M2_ind;

//cg bead radii of kinesin
vector<double> KIN_radii;

//kinesin rigid body frames
RBobj KIN_M1_RB;
RBobj KIN_M2_RB;

//kinesin rigid body forces
mvector KIN_M1_RB_force;
mvector KIN_M2_RB_force;

//kinesin rigid body torques
mvector KIN_M1_RB_torque;
mvector KIN_M2_RB_torque;

//kinesin rigid body random forces
mvector KIN_M1_RB_rand_f1;
mvector KIN_M1_RB_rand_f2;
mvector KIN_M2_RB_rand_f1;
mvector KIN_M2_RB_rand_f2;

//kinesin rigid body random torques
mvector KIN_M1_RB_rand_t1;
mvector KIN_M1_RB_rand_t2;
mvector KIN_M2_RB_rand_t1;
mvector KIN_M2_RB_rand_t2;

//kinesin diffusion tensor
double KIN_DT[6][6];

//kinesin diffusion tensor Cholesky decomposition
double KIN_DT_Chol[6][6];

//kinesin bead forces
vector<mvector> KIN_M1_force;
vector<mvector> KIN_M2_force;

//kinesin pseudo-bead forces
vector<mvector> KIN_M1_dock_cb_force;
vector<mvector> KIN_M1_dock_nl_force;
vector<mvector> KIN_M2_dock_cb_force;
vector<mvector> KIN_M2_dock_nl_force;

//kinesin bead random fluctuations
vector<mvector> KIN_M1_fluc1;
vector<mvector> KIN_M1_fluc2;
vector<mvector> KIN_M2_fluc1;
vector<mvector> KIN_M2_fluc2;

//kinesin bead force integration coefficients
vector<double> KIN_int_fc;

//kinesin bead random integration coefficients
vector<double> KIN_int_rc;

//docking interaction energy (cover bundle)
double eps_cb;

//docking interaction energy (neck linker)
double eps_nl_M1;
double eps_nl_M2;
double eps_nl;

//docking interaction stiffness
double k_att;

//TBS binding interaction energies
double eps_tbs_M1;
double eps_tbs_M2;
double eps_tbs_weak;
double eps_tbs_strong;

//CC bending stiffness
double k_bend;

//cargo rotational stiffness
double k_rot;

//MT repulsion radius
double MT_rad_rep;

//simulation parameters
//number of frames to be saved
int sim_N_frames;

//saved frames counter
int sim_frame_count = 1;

//saving period (save every number of steps)
int sim_save_period;

//steps count
int sim_step_count = 1;

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
///////////////////////INITIALIZATION/////////////////////////

void KIN_initialize_objects()
{
    //number of cg beads
    KIN_N_beads = 39;

    //read cg bead positions from file
    read_cg_beads(KIN_beads_init,"CGmodel-Kin-coord.txt");

    //update kinesin bead positions
    KIN_M1_beads = KIN_beads_init;
    KIN_M2_beads = KIN_beads_init;

    //initialize pseudo-beads
    KIN_dock_pos_cb_init.assign(6,set_vector(0,0,0));
    KIN_dock_pos_nl_init.assign(12,set_vector(0,0,0));

    for(int i=0;i<6;i++)
    {
        KIN_dock_pos_cb_init[i] = KIN_beads_init[i];
    }

    for(int i=0;i<12;i++)
    {
        KIN_dock_pos_nl_init[i] = KIN_beads_init[i+KIN_N_beads-12];
    }

    KIN_M1_dock_pos_cb = KIN_dock_pos_cb_init;
    KIN_M1_dock_pos_nl = KIN_dock_pos_nl_init;
    KIN_M2_dock_pos_cb = KIN_dock_pos_cb_init;
    KIN_M2_dock_pos_nl = KIN_dock_pos_nl_init;

    //read kinesin cg bead radii from file
    read_cg_radii(KIN_radii,"CGmodel-Kin-radii.txt");

    //set kinesin rigid body frame initial configuration
    KIN_M1_RB.pos = set_vector(111.279,14.8916,147.997);
    KIN_M1_RB.rot = set_vector(-0.633105, 1.08496, -0.983211);

    KIN_M2_RB.pos = set_vector(111.279+81.283,14.8916,147.997);
    KIN_M2_RB.rot = set_vector(-0.633105, 1.08496, -0.983211);

    //update kinesin bead coordinates
    update_rbd_coordinates(KIN_M1_beads,KIN_beads_init,KIN_M1_RB.pos,KIN_M1_RB.rot);
    update_rbd_coordinates(KIN_M2_beads,KIN_beads_init,KIN_M2_RB.pos,KIN_M2_RB.rot);

    //initialize the TBS pseudo-bead positions
    MT_TBS_M1_pos.assign(6,set_vector(0,0,0));
    MT_TBS_M2_pos.assign(6,set_vector(0,0,0));

    //read kinesin diffusion tensor from file
    read_diff_ten("DiffTen-2KIN.dat",KIN_DT);

    //read kinesin diffusion tensor Cholesky decomposition from file
    read_diff_ten("DiffTenChol-2KIN.dat",KIN_DT_Chol);

    //initialize forces and fluctuations
    KIN_M1_force.assign(KIN_N_beads,set_vector(0,0,0));
    KIN_M1_fluc1.assign(KIN_N_beads,set_vector(0,0,0));
    KIN_M1_fluc2.assign(KIN_N_beads,set_vector(0,0,0));

    KIN_M2_force.assign(KIN_N_beads,set_vector(0,0,0));
    KIN_M2_fluc1.assign(KIN_N_beads,set_vector(0,0,0));
    KIN_M2_fluc2.assign(KIN_N_beads,set_vector(0,0,0));

    KIN_M1_dock_cb_force.assign(6,set_vector(0,0,0));
    KIN_M1_dock_nl_force.assign(12,set_vector(0,0,0));
    KIN_M2_dock_cb_force.assign(6,set_vector(0,0,0));
    KIN_M2_dock_nl_force.assign(12,set_vector(0,0,0));

    //initialize integration coefficients
    KIN_int_fc.assign(KIN_N_beads,0.0);
    KIN_int_rc.assign(KIN_N_beads,0.0);

    for(int i=0;i<KIN_N_beads;i++)
    {
        KIN_int_fc[i] = h/(6*PI*eta*KIN_radii[i]);
        KIN_int_rc[i] = sqrt(2.0*h*KT/(6*PI*eta*KIN_radii[i]));
    }

    MT_rad_rep = 127.00; //A
}

//defining a function to read simulation parameters from file
void KIN_read_param()
{
    FILE * f1;
	char Ch1[100];

	f1 = fopen ("kin-param","r");

    cout << "**************************************************" << endl;
    cout << "Extracting Parameters from File" << endl;
    cout << "**************************************************" << endl;

	while (!feof(f1))
    {
        fscanf(f1,"%s\n",Ch1);

        //read the docking/binding stiffness parameter
        if (!strcmp(Ch1,"k-stiffness"))
		{
			fscanf(f1,"%lf\n",&k_att);
			cout << "Docking/Binding Stiffness: " << k_att << " A^-2" << endl;
		}

		//read the cover bundle docking energy parameter
        if (!strcmp(Ch1,"epsilon-CB"))
		{
			fscanf(f1,"%lf\n",&eps_cb);
			cout << "Cover Bundle Docking Energy: " << eps_cb << " KCal/Mol" << endl;

			eps_cb *= KCALPMOLA_TO_PN; //convert Kcal per mol to pN*A units
		}

		//read the neck linker docking energy parameter
        if (!strcmp(Ch1,"epsilon-NL"))
		{
			fscanf(f1,"%lf\n",&eps_nl);
			cout << "NL Docking Energy: " << eps_nl << " KCal/Mol" << endl;

			eps_nl *= KCALPMOLA_TO_PN; //convert Kcal per mol to pN*A units

			//default parameter assignment
			eps_nl_M1 = 0.0;
			eps_nl_M2 = eps_nl;
		}

		//read the TBS binding energy parameter
        if (!strcmp(Ch1,"epsilon-TBS-WEAK"))
		{
			fscanf(f1,"%lf\n",&eps_tbs_weak);
			cout << "Weak TBS Binding Energy: " << eps_tbs_weak << " KCal/Mol" << endl;

			eps_tbs_weak *= KCALPMOLA_TO_PN; //convert Kcal per mol to pN*A units

			//default parameter assignment
			eps_tbs_M1 = eps_tbs_weak;
		}

		//read the TBS binding energy parameter
        if (!strcmp(Ch1,"epsilon-TBS-STRONG"))
		{
			fscanf(f1,"%lf\n",&eps_tbs_strong);
			cout << "Strong TBS Binding Energy: " << eps_tbs_strong << " KCal/Mol" << endl;

			eps_tbs_strong *= KCALPMOLA_TO_PN; //convert Kcal per mol to pN*A units

			//default parameter assignment
			eps_tbs_M2 = eps_tbs_strong;
		}

        //read the CC bending stiffness parameter
        if (!strcmp(Ch1,"k-bend"))
		{
			fscanf(f1,"%lf\n",&k_bend);
			cout << "CC Bending Stiffness: " << k_bend << " KCal/Mol*rad-2" << endl;

			k_bend *= KCALPMOLA_TO_PN; //convert to pN*A
		}

		//read the cargo rotation stiffness parameter
        if (!strcmp(Ch1,"k-rotation"))
		{
			fscanf(f1,"%lf\n",&k_rot);
			cout << "Cargo Rotational Stiffness: " << k_rot << " pN*nm*rad-2" << endl;

			k_rot *= 10.0; //convert to pN*A
		}

		//read number of frames to be saved
        if (!strcmp(Ch1,"frames"))
		{
			fscanf(f1,"%d\n",&sim_N_frames);
			cout << "Number of Frames: " << sim_N_frames << " frames" << endl;
		}

		//read saving period
        if (!strcmp(Ch1,"save-period"))
		{
			fscanf(f1,"%d\n",&sim_save_period);
			cout << "Save frame every: " << sim_save_period << " steps" << endl;
		}
    }

    cout << "--------------------------------------------------" << endl << endl;
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
////////////////////////////LISTS/////////////////////////////

//motor internal repulsion interaction list
vector<index_pair> KIN_MOT_rep_list;

//motor intermolecular repulsion interaction list
vector<index_pair> KIN_M12_rep_list;

//coiled coil and motor repulsion interaction list
vector<index_pair> KIN_CC_rep_list;

//coiled coil neck linker junction and motor repulsion interaction list
vector<index_pair> KIN_CCNLJ_rep_list;

//defining a function to initiate the system interaction lists
void KIN_initialize_lists()
{
    index_pair tmp_pair;

    //initializing motor internal repulsion interaction list
    for(int i=0;i<KIN_N_beads-2;i++)
    {
        for(int j=i+2;j<KIN_N_beads;j++)
        {
            //cover bundle - self
            if((i>=0&&i<=6)&&(j>=0&&j<=6))
            {
                tmp_pair.a = i;
                tmp_pair.b = j;

                KIN_MOT_rep_list.push_back(tmp_pair);
            }

            //neck linker - self
            if((i>=26&&i<=KIN_N_beads)&&(j>=26&&j<=KIN_N_beads))
            {
                tmp_pair.a = i;
                tmp_pair.b = j;

                KIN_MOT_rep_list.push_back(tmp_pair);
            }

            //cover bundle - motor
            if((i>=0&&i<=5)&&(j>=7&&j<=KIN_N_beads))
            {
                tmp_pair.a = i;
                tmp_pair.b = j;

                KIN_MOT_rep_list.push_back(tmp_pair);
            }

            //neck linker - motor
            if((i>=0&&i<=25)&&(j>=27&&j<=KIN_N_beads))
            {
                tmp_pair.a = i;
                tmp_pair.b = j;

                KIN_MOT_rep_list.push_back(tmp_pair);
            }
        }
    }

    //initializing motor - motor repulsion interaction list
    for(int i=0;i<KIN_N_beads;i++)
    {
        for(int j=0;j<KIN_N_beads;j++)
        {
            tmp_pair.a = i;
            tmp_pair.b = j;

            if(!(i==KIN_N_beads-1&&j==KIN_N_beads-1)) KIN_M12_rep_list.push_back(tmp_pair);
        }
    }

    //initializing motor - coiled coil repulsion interaction list
    for(int i=0;i<KIN_N_beads-1;i++)
    {
        for(int j=0;j<CC_N_beads;j++)
        {
            tmp_pair.a = i;
            tmp_pair.b = j;

            KIN_CC_rep_list.push_back(tmp_pair);
        }
    }

    //initializing motor - coiled coil neck linker junction repulsion interaction list
    for(int i=0;i<KIN_N_beads-1;i++)
    {
        tmp_pair.a = i;
        tmp_pair.b = 0;

        KIN_CCNLJ_rep_list.push_back(tmp_pair);
    }
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
///////////////////////NEIGHBOR LISTS/////////////////////////

//motor internal repulsion neighbor list
vector<index_pair> KIN_M1_rep_Nlist;
vector<index_pair> KIN_M2_rep_Nlist;

//motor inter motor repulsion neighbor list
vector<index_pair> KIN_M12_rep_Nlist;

//motor - coiled coil repulsion neighbor list
vector<index_pair> KIN_M1_CC_rep_Nlist;
vector<index_pair> KIN_M2_CC_rep_Nlist;

//motor - coiled coil neck linker junction repulsion neighbor list
vector<index_pair> KIN_M1_CCNLJ_rep_Nlist;
vector<index_pair> KIN_M2_CCNLJ_rep_Nlist;

//position records
vector<mvector> KIN_dr_pos0,KIN_dr_pos1;

//displacement record
vector<double> KIN_dr;

//defining a function to create the position and displacement records
void KIN_create_pos_record()
{
    int In1;

    //summing the total number of particles involved in neighbor lists
    In1 = 2*KIN_N_beads;

    //creating position records
    KIN_dr_pos0.assign(In1,set_vector(0,0,0));
    KIN_dr_pos1.assign(In1,set_vector(0,0,0));

    //creating displacement record
    KIN_dr.assign(In1,0.0);
}

//defining a function to reset the displacement record
void KIN_reset_displacement_record()
{
    int In1;

    //summing the total number of particles involved in neighbor lists
    In1 = 2*KIN_N_beads;

    for(int i=0;i<In1;i++)
    {
        KIN_dr[i] = 0.0;
    }
}

//defining a function to update the displacement record
void KIN_update_displacement_record()
{
    int In1;
    mvector dr;

    //summing the total number of particles involved in neighbor lists
    In1 = 2*KIN_N_beads;

    for(int i=0;i<In1;i++)
    {
        dr = KIN_dr_pos1[i] - KIN_dr_pos0[i];
        KIN_dr[i] = dr.square();
    }
}

//defining a function to reset the position record
void KIN_reset_pos_record()
{
    int In1 = 0;

    //M1 beads
    for(int i=0;i<KIN_N_beads;i++)
    {
        KIN_dr_pos0[i+In1] = KIN_M1_beads[i];
    }

    In1 += KIN_N_beads;

    //M2 beads
    for(int i=0;i<KIN_N_beads;i++)
    {
        KIN_dr_pos0[i+In1] = KIN_M2_beads[i];
    }

    In1 += KIN_N_beads;

    //reset record
    for(int i=0;i<In1;i++)
    {
        KIN_dr_pos1[i] = KIN_dr_pos0[i];
    }
}

//defining a function to update the position record
void KIN_update_pos_record()
{
    int In1 = 0;

    //M1 beads
    for(int i=0;i<KIN_N_beads;i++)
    {
        KIN_dr_pos1[i+In1] = KIN_M1_beads[i];
    }

    In1 += KIN_N_beads;

    //M2 beads
    for(int i=0;i<KIN_N_beads;i++)
    {
        KIN_dr_pos1[i+In1] = KIN_M2_beads[i];
    }

    In1 += KIN_N_beads;
}

//defining a function to update the kinesin neighbor lists
void KIN_update_neighbor_lists()
{
    //clearing lists
    KIN_M1_rep_Nlist.clear();
    KIN_M2_rep_Nlist.clear();
    KIN_M12_rep_Nlist.clear();
    KIN_M1_CC_rep_Nlist.clear();
    KIN_M2_CC_rep_Nlist.clear();
    KIN_M1_CCNLJ_rep_Nlist.clear();
    KIN_M2_CCNLJ_rep_Nlist.clear();

    //updating lists
    KIN_M1_rep_Nlist = create_neighbor_list(KIN_M1_beads,KIN_MOT_rep_list);
    KIN_M2_rep_Nlist = create_neighbor_list(KIN_M2_beads,KIN_MOT_rep_list);
    KIN_M12_rep_Nlist = create_neighbor_list(KIN_M1_beads,KIN_M2_beads,KIN_M12_rep_list);
    KIN_M1_CC_rep_Nlist = create_neighbor_list(KIN_M1_beads,CC_beads,KIN_CC_rep_list);
    KIN_M2_CC_rep_Nlist = create_neighbor_list(KIN_M2_beads,CC_beads,KIN_CC_rep_list);
    KIN_M1_CCNLJ_rep_Nlist = create_neighbor_list(KIN_M1_beads,CCNLJ_beads,KIN_CCNLJ_rep_list);
    KIN_M2_CCNLJ_rep_Nlist = create_neighbor_list(KIN_M2_beads,CCNLJ_beads,KIN_CCNLJ_rep_list);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
///////////////////////BAOAB DYNAMICS/////////////////////////

//defining a function to initialize the random forces in the system
void KIN_initialize_BAOAB()
{
    //get rigid body random forces
    KIN_M1_RB_rand_f1 = get_rand_mvector(init);
    KIN_M1_RB_rand_f2 = get_rand_mvector(init);
    KIN_M2_RB_rand_f1 = get_rand_mvector(init);
    KIN_M2_RB_rand_f2 = get_rand_mvector(init);
    CCNLJ_RB_rand_f1 = get_rand_mvector(init);
    CCNLJ_RB_rand_f2 = get_rand_mvector(init);
    CAR_RB_rand_f1 = get_rand_mvector(init);
    CAR_RB_rand_f2 = get_rand_mvector(init);

    //get rigid body random torques
    KIN_M1_RB_rand_t1 = get_rand_mvector(init);
    KIN_M1_RB_rand_t2 = get_rand_mvector(init);
    KIN_M2_RB_rand_t1 = get_rand_mvector(init);
    KIN_M2_RB_rand_t2 = get_rand_mvector(init);
    CCNLJ_RB_rand_t1 = get_rand_mvector(init);
    CCNLJ_RB_rand_t2 = get_rand_mvector(init);
    CAR_RB_rand_t1 = get_rand_mvector(init);
    CAR_RB_rand_t2 = get_rand_mvector(init);

    //correlate forces and torques
    correlate_force_torque(KIN_M1_RB_rand_f1,KIN_M1_RB_rand_t1,KIN_M1_RB.rot,KIN_DT_Chol);
    correlate_force_torque(KIN_M1_RB_rand_f2,KIN_M1_RB_rand_t2,KIN_M1_RB.rot,KIN_DT_Chol);
    correlate_force_torque(KIN_M2_RB_rand_f1,KIN_M2_RB_rand_t1,KIN_M2_RB.rot,KIN_DT_Chol);
    correlate_force_torque(KIN_M2_RB_rand_f2,KIN_M2_RB_rand_t2,KIN_M2_RB.rot,KIN_DT_Chol);
    correlate_force_torque(CCNLJ_RB_rand_f1,CCNLJ_RB_rand_t1,CCNLJ_RB.rot,CCNLJ_DT_Chol);
    correlate_force_torque(CCNLJ_RB_rand_f2,CCNLJ_RB_rand_t2,CCNLJ_RB.rot,CCNLJ_DT_Chol);
    correlate_force_torque(CAR_RB_rand_f1,CAR_RB_rand_t1,CAR_RB.rot,CAR_DT_Chol);
    correlate_force_torque(CAR_RB_rand_f2,CAR_RB_rand_t2,CAR_RB.rot,CAR_DT_Chol);

    //get random forces
    for(int i=0;i<6;i++)
    {
        KIN_M1_fluc1[i] = get_rand_mvector(init);
        KIN_M1_fluc2[i] = get_rand_mvector(init);

        KIN_M2_fluc1[i] = get_rand_mvector(init);
        KIN_M2_fluc2[i] = get_rand_mvector(init);
    }

    for(int i=27;i<KIN_N_beads;i++)
    {
        KIN_M1_fluc1[i] = get_rand_mvector(init);
        KIN_M1_fluc2[i] = get_rand_mvector(init);

        KIN_M2_fluc1[i] = get_rand_mvector(init);
        KIN_M2_fluc2[i] = get_rand_mvector(init);
    }

    for(int i=0;i<CC_N_beads;i++)
    {
        CC_fluc1[i] = get_rand_mvector(init);
        CC_fluc2[i] = get_rand_mvector(init);
    }
}

void KIN_BAOAB_get_fluc()
{
    //keep rigid body random forces from previous step
    KIN_M1_RB_rand_f1 = KIN_M1_RB_rand_f2;
    KIN_M2_RB_rand_f1 = KIN_M2_RB_rand_f2;
    CCNLJ_RB_rand_f1 = CCNLJ_RB_rand_f2;
    CAR_RB_rand_f1 = CAR_RB_rand_f2;

    //keep rigid body random torques from previous step
    KIN_M1_RB_rand_t1 = KIN_M1_RB_rand_t2;
    KIN_M2_RB_rand_t1 = KIN_M2_RB_rand_t2;
    CCNLJ_RB_rand_t1 = CCNLJ_RB_rand_t2;
    CAR_RB_rand_t1 = CAR_RB_rand_t2;

    //get rigid body random forces
    KIN_M1_RB_rand_f2 = get_rand_mvector(init);
    KIN_M2_RB_rand_f2 = get_rand_mvector(init);
    CCNLJ_RB_rand_f2 = get_rand_mvector(init);
    CAR_RB_rand_f2 = get_rand_mvector(init);

    //get rigid body random torques
    KIN_M1_RB_rand_t2 = get_rand_mvector(init);
    KIN_M2_RB_rand_t2 = get_rand_mvector(init);
    CCNLJ_RB_rand_t2 = get_rand_mvector(init);
    CAR_RB_rand_t2 = get_rand_mvector(init);

    //correlate forces and torques
    correlate_force_torque(KIN_M1_RB_rand_f2,KIN_M1_RB_rand_t2,KIN_M1_RB.rot,KIN_DT_Chol);
    correlate_force_torque(KIN_M2_RB_rand_f2,KIN_M2_RB_rand_t2,KIN_M2_RB.rot,KIN_DT_Chol);
    correlate_force_torque(CCNLJ_RB_rand_f2,CCNLJ_RB_rand_t2,CCNLJ_RB.rot,CCNLJ_DT_Chol);
    correlate_force_torque(CAR_RB_rand_f2,CAR_RB_rand_t2,CAR_RB.rot,CAR_DT_Chol);

    //keep random forces from previous step
    KIN_M1_fluc1 = KIN_M1_fluc2;
    KIN_M2_fluc1 = KIN_M2_fluc2;
    CC_fluc1 = CC_fluc2;

    //get random forces
    for(int i=0;i<6;i++)
    {
        KIN_M1_fluc2[i] = get_rand_mvector(init);
        KIN_M2_fluc2[i] = get_rand_mvector(init);
    }

    for(int i=27;i<KIN_N_beads;i++)
    {
        KIN_M1_fluc2[i] = get_rand_mvector(init);
        KIN_M2_fluc2[i] = get_rand_mvector(init);
    }

    for(int i=0;i<CC_N_beads;i++)
    {
        CC_fluc2[i] = get_rand_mvector(init);
    }
}

//defining a function to integrate the time step for all components in the system
void KIN_BAOAB_integrate()
{
    //integrate kinesin rigid bodies
    RBD_integrate_brownian(KIN_M1_RB,KIN_M1_RB_force,KIN_M1_RB_torque,0.5*(KIN_M1_RB_rand_f1+KIN_M1_RB_rand_f2),0.5*(KIN_M1_RB_rand_t1+KIN_M1_RB_rand_t2));
    RBD_integrate_brownian(KIN_M2_RB,KIN_M2_RB_force,KIN_M2_RB_torque,0.5*(KIN_M2_RB_rand_f1+KIN_M2_RB_rand_f2),0.5*(KIN_M2_RB_rand_t1+KIN_M2_RB_rand_t2));
    RBD_integrate_brownian(CCNLJ_RB,CCNLJ_RB_force,CCNLJ_RB_torque,0.5*(CCNLJ_RB_rand_f1+CCNLJ_RB_rand_f2),0.5*(CCNLJ_RB_rand_t1+CCNLJ_RB_rand_t2));
    RBD_integrate_brownian(CAR_RB,CAR_RB_force,CAR_RB_torque,0.5*(CAR_RB_rand_f1+CAR_RB_rand_f2),0.5*(CAR_RB_rand_t1+CAR_RB_rand_t2));

    //update bead coordinates of rigid bodies
    update_rbd_coordinates(KIN_M1_beads,KIN_beads_init,KIN_M1_RB.pos,KIN_M1_RB.rot,6,27);
    update_rbd_coordinates(KIN_M2_beads,KIN_beads_init,KIN_M2_RB.pos,KIN_M2_RB.rot,6,27);
    update_rbd_coordinates(CCNLJ_beads,CCNLJ_beads_init,CCNLJ_RB.pos,CCNLJ_RB.rot);
    update_rbd_coordinates(CAR_beads,CAR_beads_init,CAR_RB.pos,CAR_RB.rot);

    //update pseudo-bead coordinates of rigid bodies
    update_rbd_coordinates(KIN_M1_dock_pos_cb,KIN_dock_pos_cb_init,KIN_M1_RB.pos,KIN_M1_RB.rot,0,6);
    update_rbd_coordinates(KIN_M1_dock_pos_nl,KIN_dock_pos_nl_init,KIN_M1_RB.pos,KIN_M1_RB.rot,0,12);
    update_rbd_coordinates(KIN_M2_dock_pos_cb,KIN_dock_pos_cb_init,KIN_M2_RB.pos,KIN_M2_RB.rot,0,6);
    update_rbd_coordinates(KIN_M2_dock_pos_nl,KIN_dock_pos_nl_init,KIN_M2_RB.pos,KIN_M2_RB.rot,0,12);

    KIN_M1_beads[KIN_N_beads-1] = CCNLJ_beads[1];
    KIN_M2_beads[KIN_N_beads-1] = CCNLJ_beads[2];

    //integrate beads
    integrate_brownian_BAOAB(KIN_M1_beads,KIN_M1_force,KIN_M1_fluc1,KIN_M1_fluc2,KIN_int_fc,KIN_int_rc,0,6);
    integrate_brownian_BAOAB(KIN_M2_beads,KIN_M2_force,KIN_M2_fluc1,KIN_M2_fluc2,KIN_int_fc,KIN_int_rc,0,6);

    integrate_brownian_BAOAB(KIN_M1_beads,KIN_M1_force,KIN_M1_fluc1,KIN_M1_fluc2,KIN_int_fc,KIN_int_rc,27,KIN_N_beads-1);
    integrate_brownian_BAOAB(KIN_M2_beads,KIN_M2_force,KIN_M2_fluc1,KIN_M2_fluc2,KIN_int_fc,KIN_int_rc,27,KIN_N_beads-1);

    integrate_brownian_BAOAB(CC_beads,CC_force,CC_fluc1,CC_fluc2,CC_int_fc,CC_int_rc,0,CC_N_beads);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////MECHANICS///////////////////////////

//defining a function to resolve the forces and torques in the system
void KIN_resolve_forces_torques()
{
    //obtain resultant forces
    KIN_M1_RB_force += calc_RB_resultant_force(KIN_M1_force,6,27);
    KIN_M2_RB_force += calc_RB_resultant_force(KIN_M2_force,6,27);
    CCNLJ_RB_force += calc_RB_resultant_force(CCNLJ_force,0,4);
    CAR_RB_force += calc_RB_resultant_force(CAR_force,0,2);

    KIN_M1_RB_force += calc_RB_resultant_force(KIN_M1_dock_cb_force,0,6);
    KIN_M1_RB_force += calc_RB_resultant_force(KIN_M1_dock_nl_force,0,12);
    KIN_M2_RB_force += calc_RB_resultant_force(KIN_M2_dock_cb_force,0,6);
    KIN_M2_RB_force += calc_RB_resultant_force(KIN_M2_dock_nl_force,0,12);

    //obtain resultant torques
    KIN_M1_RB_torque += calc_RB_resultant_torque(KIN_M1_force,KIN_M1_beads,KIN_M1_RB.pos,6,27);
    KIN_M2_RB_torque += calc_RB_resultant_torque(KIN_M2_force,KIN_M2_beads,KIN_M2_RB.pos,6,27);
    CCNLJ_RB_torque += calc_RB_resultant_torque(CCNLJ_force,CCNLJ_beads,CCNLJ_RB.pos,0,4);
    CAR_RB_torque += calc_RB_resultant_torque(CAR_force,CAR_beads,CAR_RB.pos,0,2);

    KIN_M1_RB_torque += calc_RB_resultant_torque(KIN_M1_dock_cb_force,KIN_M1_dock_pos_cb,KIN_M1_RB.pos,0,6);
    KIN_M1_RB_torque += calc_RB_resultant_torque(KIN_M1_dock_nl_force,KIN_M1_dock_pos_nl,KIN_M1_RB.pos,0,12);
    KIN_M2_RB_torque += calc_RB_resultant_torque(KIN_M2_dock_cb_force,KIN_M2_dock_pos_cb,KIN_M2_RB.pos,0,6);
    KIN_M2_RB_torque += calc_RB_resultant_torque(KIN_M2_dock_nl_force,KIN_M2_dock_pos_nl,KIN_M2_RB.pos,0,12);

    //correlate forces and torques
    correlate_force_torque(KIN_M1_RB_force,KIN_M1_RB_torque,KIN_M1_RB.rot,KIN_DT);
    correlate_force_torque(KIN_M2_RB_force,KIN_M2_RB_torque,KIN_M2_RB.rot,KIN_DT);
    correlate_force_torque(CCNLJ_RB_force,CCNLJ_RB_torque,CCNLJ_RB.rot,CCNLJ_DT);
    correlate_force_torque(CAR_RB_force,CAR_RB_torque,CAR_RB.rot,CAR_DT);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
///////////////////////////FORCES/////////////////////////////

//defining a function to apply the connectivity forces
void KIN_apply_forces_connectivity()
{
    mvector tmp_f;

    //cover bundle connectivity forces
    for(int i=0;i<6;i++)
    {
        //M1
        tmp_f = calc_fene(KIN_M1_beads[i],KIN_M1_beads[i+1],3.8);

        KIN_M1_force[i] += tmp_f;
        KIN_M1_force[i+1] -= tmp_f;

        //M2
        tmp_f = calc_fene(KIN_M2_beads[i],KIN_M2_beads[i+1],3.8);

        KIN_M2_force[i] += tmp_f;
        KIN_M2_force[i+1] -= tmp_f;
    }

    //neck linker connectivity forces
    for(int i=26;i<KIN_N_beads-1;i++)
    {
        //M1
        tmp_f = calc_fene(KIN_M1_beads[i],KIN_M1_beads[i+1],3.8);

        KIN_M1_force[i] += tmp_f;
        KIN_M1_force[i+1] -= tmp_f;

        //M2
        tmp_f = calc_fene(KIN_M2_beads[i],KIN_M2_beads[i+1],3.8);

        KIN_M2_force[i] += tmp_f;
        KIN_M2_force[i+1] -= tmp_f;
    }

    //coiled coil connectivity forces
    for(int i=0;i<CC_N_beads-1;i++)
    {
        tmp_f = calc_fene(CC_beads[i],CC_beads[i+1],2.0*CC_radii);

        CC_force[i] += tmp_f;
        CC_force[i+1] -= tmp_f;
    }

    //CC-NL-Junction
    tmp_f = calc_fene(CCNLJ_beads[0],CC_beads[0],2.0*CC_radii);

    CCNLJ_force[0] += tmp_f;
    CC_force[0] -= tmp_f;
}

//defining a function to apply the excluded volume forces within molecules
void KIN_apply_forces_excluded_volume_internal()
{
    mvector tmp_f;
    int In1,ind1,ind2;

    //M1
    In1 = (int) KIN_M1_rep_Nlist.size();

    for(int i=0;i<In1;i++)
    {
        ind1 = KIN_M1_rep_Nlist[i].a;
        ind2 = KIN_M1_rep_Nlist[i].b;

        //M1
        tmp_f = calc_repulsion(KIN_M1_beads[ind1],KIN_M1_beads[ind2],KIN_radii[ind1],KIN_radii[ind2],3.8);

        KIN_M1_force[ind1] += tmp_f;
        KIN_M1_force[ind2] -= tmp_f;
    }

    //M2
    In1 = (int) KIN_M2_rep_Nlist.size();

    for(int i=0;i<In1;i++)
    {
        ind1 = KIN_M2_rep_Nlist[i].a;
        ind2 = KIN_M2_rep_Nlist[i].b;

        //M2
        tmp_f = calc_repulsion(KIN_M2_beads[ind1],KIN_M2_beads[ind2],KIN_radii[ind1],KIN_radii[ind2],3.8);

        KIN_M2_force[ind1] += tmp_f;
        KIN_M2_force[ind2] -= tmp_f;
    }
}

//defining a function to apply the excluded volume forces between molecules
void KIN_apply_forces_excluded_volume_external()
{
    mvector tmp_f;
    int In1,ind1,ind2;

    In1 = (int) KIN_M12_rep_Nlist.size();

    //motor external excluded volume forces
    for(int i=0;i<In1;i++)
    {
        ind1 = KIN_M12_rep_Nlist[i].a;
        ind2 = KIN_M12_rep_Nlist[i].b;

        tmp_f = calc_repulsion(KIN_M1_beads[ind1],KIN_M2_beads[ind2],KIN_radii[ind1],KIN_radii[ind2],3.8);

        KIN_M1_force[ind1] += tmp_f;
        KIN_M2_force[ind2] -= tmp_f;
    }

    In1 = (int) KIN_M1_CC_rep_Nlist.size();

    //M1-CC external excluded volume forces
    for(int i=0;i<In1;i++)
    {
        ind1 = KIN_M1_CC_rep_Nlist[i].a;
        ind2 = KIN_M1_CC_rep_Nlist[i].b;

        tmp_f = calc_repulsion(KIN_M1_beads[ind1],CC_beads[ind2],KIN_radii[ind1],CC_radii,3.8);

        KIN_M1_force[ind1] += tmp_f;
        CC_force[ind2] -= tmp_f;
    }

    In1 = (int) KIN_M2_CC_rep_Nlist.size();

    //M1-CC external excluded volume forces
    for(int i=0;i<In1;i++)
    {
        ind1 = KIN_M2_CC_rep_Nlist[i].a;
        ind2 = KIN_M2_CC_rep_Nlist[i].b;

        tmp_f = calc_repulsion(KIN_M2_beads[ind1],CC_beads[ind2],KIN_radii[ind1],CC_radii,3.8);

        KIN_M2_force[ind1] += tmp_f;
        CC_force[ind2] -= tmp_f;
    }

    In1 = (int) KIN_M1_CCNLJ_rep_Nlist.size();

    //M1-CCNLJ external excluded volume forces
    for(int i=0;i<In1;i++)
    {
        ind1 = KIN_M1_CCNLJ_rep_Nlist[i].a;
        ind2 = KIN_M1_CCNLJ_rep_Nlist[i].b;

        tmp_f = calc_repulsion(KIN_M1_beads[ind1],CCNLJ_beads[ind2],KIN_radii[ind1],CC_radii,3.8);

        KIN_M1_force[ind1] += tmp_f;
        CCNLJ_force[ind2] -= tmp_f;
    }

    In1 = (int) KIN_M2_CCNLJ_rep_Nlist.size();

    //M2-CCNLJ external excluded volume forces
    for(int i=0;i<In1;i++)
    {
        ind1 = KIN_M2_CCNLJ_rep_Nlist[i].a;
        ind2 = KIN_M2_CCNLJ_rep_Nlist[i].b;

        tmp_f = calc_repulsion(KIN_M2_beads[ind1],CCNLJ_beads[ind2],KIN_radii[ind1],CC_radii,3.8);

        KIN_M2_force[ind1] += tmp_f;
        CCNLJ_force[ind2] -= tmp_f;
    }
}

//defining a function to apply angular forces
void CC_apply_forces_angular()
{
    mvector tmp_f1,tmp_f2,Mv1,Mv2;

    //angular forces along the coiled coil
    for(int i=0;i<CC_N_beads-2;i++)
    {
        //defining the mid points
        Mv1 = CC_beads[i] - CC_beads[i+1];
        Mv2 = CC_beads[i+2] - CC_beads[i+1];

        calc_angular(Mv1,Mv2,tmp_f1,tmp_f2,PI,k_bend);

        CC_force[i] += tmp_f1;

        CC_force[i+1] -= tmp_f1 + tmp_f2;

        CC_force[i+2] += tmp_f2;
    }

    //CC-NL-Junction
    //defining the mid points
    Mv1 = CCNLJ_beads[0] - CC_beads[0];
    Mv2 = CC_beads[1] - CC_beads[0];

    calc_angular(Mv1,Mv2,tmp_f1,tmp_f2,PI,k_bend);

    CCNLJ_force[0] += tmp_f1;

    CC_force[0] -= tmp_f1 + tmp_f2;

    CC_force[1] += tmp_f2;

    //defining the mid points
    Mv1 = CCNLJ_beads[3] - CCNLJ_beads[0];
    Mv2 = CC_beads[0] - CCNLJ_beads[0];

    calc_angular(Mv1,Mv2,tmp_f1,tmp_f2,0.0,k_bend);

    CCNLJ_force[3] += tmp_f1;

    CCNLJ_force[0] -= tmp_f1 + tmp_f2;

    CC_force[0] += tmp_f2;
}

//defining a function to apply docking forces
void KIN_apply_forces_docking()
{
    mvector tmp_f;

    //cover bundle docking forces
    for(int i=0;i<6;i++)
    {
        //M1
        tmp_f = calc_attraction(KIN_M1_beads[i],KIN_M1_dock_pos_cb[i],k_att,eps_cb);

        KIN_M1_force[i] += tmp_f;
        KIN_M1_dock_cb_force[i] -= tmp_f;

        //M2
        tmp_f = calc_attraction(KIN_M2_beads[i],KIN_M2_dock_pos_cb[i],k_att,eps_cb);

        KIN_M2_force[i] += tmp_f;
        KIN_M2_dock_cb_force[i] -= tmp_f;
    }

    //neck linker docking forces
    for(int i=0;i<12;i++)
    {
        //M1
        tmp_f = calc_attraction(KIN_M1_beads[i+KIN_N_beads-12],KIN_M1_dock_pos_nl[i],k_att,eps_nl_M1);

        KIN_M1_force[i+KIN_N_beads-12] += tmp_f;
        KIN_M1_dock_nl_force[i] -= tmp_f;

        //M2
        tmp_f = calc_attraction(KIN_M2_beads[i+KIN_N_beads-12],KIN_M2_dock_pos_nl[i],k_att,eps_nl_M2);

        KIN_M2_force[i+KIN_N_beads-12] += tmp_f;
        KIN_M2_dock_nl_force[i] -= tmp_f;
    }

}

//defining a function to apply the forces involved in the neck linker - coiled coil junction
void KIN_apply_forces_CCNLJ()
{
    //convert NL ends forces to pseudo bead forces
    CCNLJ_force[1] = KIN_M1_force[KIN_N_beads-1];
    CCNLJ_force[2] = KIN_M2_force[KIN_N_beads-1];
}

//defining a function to apply the forces involved in motor - TBS binding (M1)
void MT_apply_forces_TBS_M1()
{
    mvector tmp_f[6];

    tmp_f[0] = calc_attraction(KIN_M1_beads[18],MT_TBS_M1_pos[0],k_att,eps_tbs_M1/6.0);
    tmp_f[1] = calc_attraction(KIN_M1_beads[19],MT_TBS_M1_pos[1],k_att,eps_tbs_M1/6.0);
    tmp_f[2] = calc_attraction(KIN_M1_beads[20],MT_TBS_M1_pos[2],k_att,eps_tbs_M1/6.0);
    tmp_f[3] = calc_attraction(KIN_M1_beads[22],MT_TBS_M1_pos[3],k_att,eps_tbs_M1/6.0);
    tmp_f[4] = calc_attraction(KIN_M1_beads[23],MT_TBS_M1_pos[4],k_att,eps_tbs_M1/6.0);
    tmp_f[5] = calc_attraction(KIN_M1_beads[24],MT_TBS_M1_pos[5],k_att,eps_tbs_M1/6.0);

    KIN_M1_force[18] += tmp_f[0];
    KIN_M1_force[19] += tmp_f[1];
    KIN_M1_force[20] += tmp_f[2];
    KIN_M1_force[22] += tmp_f[3];
    KIN_M1_force[23] += tmp_f[4];
    KIN_M1_force[24] += tmp_f[5];
}

//defining a function to apply the forces involved in motor - TBS binding (M2)
void MT_apply_forces_TBS_M2()
{
    mvector tmp_f[6];

    tmp_f[0] = calc_attraction(KIN_M2_beads[18],MT_TBS_M2_pos[0],k_att,eps_tbs_M2/6.0);
    tmp_f[1] = calc_attraction(KIN_M2_beads[19],MT_TBS_M2_pos[1],k_att,eps_tbs_M2/6.0);
    tmp_f[2] = calc_attraction(KIN_M2_beads[20],MT_TBS_M2_pos[2],k_att,eps_tbs_M2/6.0);
    tmp_f[3] = calc_attraction(KIN_M2_beads[22],MT_TBS_M2_pos[3],k_att,eps_tbs_M2/6.0);
    tmp_f[4] = calc_attraction(KIN_M2_beads[23],MT_TBS_M2_pos[4],k_att,eps_tbs_M2/6.0);
    tmp_f[5] = calc_attraction(KIN_M2_beads[24],MT_TBS_M2_pos[5],k_att,eps_tbs_M2/6.0);

    KIN_M2_force[18] += tmp_f[0];
    KIN_M2_force[19] += tmp_f[1];
    KIN_M2_force[20] += tmp_f[2];
    KIN_M2_force[22] += tmp_f[3];
    KIN_M2_force[23] += tmp_f[4];
    KIN_M2_force[24] += tmp_f[5];
}

//defining a function to apply repulsion forces between the MT and the rest of the particles in the system
void MT_apply_forces_excluded_volume()
{
    mvector tmp_f;

    //motor
    for(int i=0;i<KIN_N_beads;i++)
    {
        //M1
        tmp_f = calc_repulsion_MT(KIN_M1_beads[i],KIN_radii[i],MT_rad_rep,3.8);

        KIN_M1_force[i] += tmp_f;

        //M2
        tmp_f = calc_repulsion_MT(KIN_M2_beads[i],KIN_radii[i],MT_rad_rep,3.8);

        KIN_M2_force[i] += tmp_f;
    }

    //CC
    for(int i=0;i<CC_N_beads;i++)
    {
        tmp_f = calc_repulsion_MT(CC_beads[i],CC_radii,MT_rad_rep,3.8);

        CC_force[i] += tmp_f;
    }

}

//defining a function to calculate the observed angle between cargo and CC-NL-Junction
double CAR_get_obs_ang()
{
    double ang1,ang2,Db1,Db2;
    mvector Mv1;
    mvector VirtAxes1_x,VirtAxes1_z;
    mvector VirtAxes2_x,VirtAxes2_z;

    VirtAxes1_x = rotate_vector_rodrigues(set_vector(1,0,0),CAR_RB.rot);
    VirtAxes2_x = rotate_vector_rodrigues(set_vector(1,0,0),CCNLJ_RB.rot);

    VirtAxes1_z = rotate_vector_rodrigues(set_vector(0,0,1),CAR_RB.rot);
    VirtAxes2_z = rotate_vector_rodrigues(set_vector(0,0,1),CCNLJ_RB.rot);

    Db1 = VirtAxes1_z*VirtAxes2_z;

    if(abs(Db1)>1.0) Db1 /= abs(Db1);

    ang1 = acos(Db1);

    Mv1 = VirtAxes2_z^VirtAxes1_z;

    if(Mv1.length()!=0) Mv1 = ang1*unit_vector(Mv1);

    VirtAxes2_x = rotate_vector_rodrigues(VirtAxes2_x,Mv1);

    Db1 = VirtAxes1_x*VirtAxes2_x;

    Db2 = (VirtAxes1_x^VirtAxes2_x)*VirtAxes1_z;

    if(abs(Db1)>1.0) Db1 /= abs(Db1);

    ang2 = acos(Db1);

    if(Db2<0) ang2 *= -1.0;

    return ang2;
}

//defining a function to apply cargo forces
void CAR_apply_forces()
{
    double ang_tmp;
    mvector tmp_f1,tmp_f2,tmp_f3,tmp_f4,Mv1,Mv2;

    //applying connectivity force between cargo and coiled coil
    tmp_f1 = calc_fene(CAR_beads[1],CC_beads.back(),2.0*CC_radii);

    CAR_force[1] += tmp_f1;
    CC_force.back() -= tmp_f1;

    //applying angular bending forces between cargo and coiled coil
    Mv1 = CC_beads[CC_N_beads-2] - CC_beads[CC_N_beads-1];
    Mv2 = CAR_beads[1] - CC_beads[CC_N_beads-1];

    calc_angular(Mv1,Mv2,tmp_f1,tmp_f2,PI,88.0*KCALPMOLA_TO_PN);

    CC_force[CC_N_beads-2] += tmp_f1;

    CC_force[CC_N_beads-1] -= tmp_f1 + tmp_f2;

    CAR_force[1] += tmp_f2;

    Mv1 = CC_beads[CC_N_beads-1] - CAR_beads[1];
    Mv2 = CAR_beads[0] - CAR_beads[1];

    calc_angular(Mv1,Mv2,tmp_f1,tmp_f2,PI,88.0*KCALPMOLA_TO_PN);

    CC_force[CC_N_beads-1] += tmp_f1;

    CAR_force[1] -= tmp_f1 + tmp_f2;

    CAR_force[0] += tmp_f2;

    //applying rotational constraint force between cargo and CC-NL-Junction
    ang_tmp = CAR_get_obs_ang();

    while(ang_tmp+CAR_obs_ang_shift-CAR_obs_ang>PI)
    {
        CAR_obs_ang_shift -= 2.0*PI;
    }

    while(ang_tmp+CAR_obs_ang_shift-CAR_obs_ang<-PI)
    {
        CAR_obs_ang_shift += 2.0*PI;
    }

    ang_tmp += CAR_obs_ang_shift;

    CAR_obs_ang = ang_tmp;

    Mv1 = rotate_vector_rodrigues(set_vector(0,0,1),CAR_RB.rot);
    Mv2 = rotate_vector_rodrigues(set_vector(0,0,1),CCNLJ_RB.rot);

    CAR_RB_torque += k_rot*CAR_obs_ang*Mv1;
    CCNLJ_RB_torque -= k_rot*CAR_obs_ang*Mv2;

    //apply repulsion force between cargo and MT
    tmp_f1 = calc_repulsion_MT(CAR_beads[0],CAR_rad,MT_rad_rep,3.8);

    CAR_force[0] += tmp_f1;

    //apply repulsion forces between cargo and the motors
    tmp_f1 = calc_repulsion(KIN_M1_RB.pos,CAR_beads[0],30.0,CAR_rad,3.8);
    tmp_f2 = calc_repulsion(KIN_M2_RB.pos,CAR_beads[0],30.0,CAR_rad,3.8);

    KIN_M1_RB_force += tmp_f1;
    CAR_force[0] -= tmp_f1;

    KIN_M2_RB_force += tmp_f2;
    CAR_force[0] -= tmp_f2;
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
///////////////////////////IO DATA////////////////////////////

//defining a function to write a snap shot of the system to file
void write_snapshot(string path)
{
    FILE * f1;

    f1 = fopen(path.c_str(),"w");

    //write position and rotation of rigid objects
    //M1
    fprintf(f1,"%lf %lf %lf\n",KIN_M1_RB.pos.x,KIN_M1_RB.pos.y,KIN_M1_RB.pos.z);
    fprintf(f1,"%lf %lf %lf\n",KIN_M1_RB.rot.x,KIN_M1_RB.rot.y,KIN_M1_RB.rot.z);

    //M2
    fprintf(f1,"%lf %lf %lf\n",KIN_M2_RB.pos.x,KIN_M2_RB.pos.y,KIN_M2_RB.pos.z);
    fprintf(f1,"%lf %lf %lf\n",KIN_M2_RB.rot.x,KIN_M2_RB.rot.y,KIN_M2_RB.rot.z);

    //CC-NL-Junction
    fprintf(f1,"%lf %lf %lf\n",CCNLJ_RB.pos.x,CCNLJ_RB.pos.y,CCNLJ_RB.pos.z);
    fprintf(f1,"%lf %lf %lf\n",CCNLJ_RB.rot.x,CCNLJ_RB.rot.y,CCNLJ_RB.rot.z);

    //cargo
    fprintf(f1,"%lf %lf %lf\n",CAR_RB.pos.x,CAR_RB.pos.y,CAR_RB.pos.z);
    fprintf(f1,"%lf %lf %lf\n",CAR_RB.rot.x,CAR_RB.rot.y,CAR_RB.rot.z);

    //write bead positions
    //M1
    for(int i=0;i<KIN_N_beads;i++)
    {
        fprintf(f1,"%lf %lf %lf\n",KIN_M1_beads[i].x,KIN_M1_beads[i].y,KIN_M1_beads[i].z);
    }

    //M2
    for(int i=0;i<KIN_N_beads;i++)
    {
        fprintf(f1,"%lf %lf %lf\n",KIN_M2_beads[i].x,KIN_M2_beads[i].y,KIN_M2_beads[i].z);
    }

    //CC
    for(int i=0;i<CC_N_beads;i++)
    {
        fprintf(f1,"%lf %lf %lf\n",CC_beads[i].x,CC_beads[i].y,CC_beads[i].z);
    }

    //CC-NL-Junction
    for(int i=0;i<4;i++)
    {
        fprintf(f1,"%lf %lf %lf\n",CCNLJ_beads[i].x,CCNLJ_beads[i].y,CCNLJ_beads[i].z);
    }

    //write bead positions of the cargo
    for(int i=0;i<2;i++)
    {
        fprintf(f1,"%lf %lf %lf\n",CAR_beads[i].x,CAR_beads[i].y,CAR_beads[i].z);
    }

    //write angular data
    fprintf(f1,"%lf %lf\n",CAR_obs_ang,CAR_obs_ang_shift);

    fclose(f1);
}

//defining a function to read a snap shot of the system from file
void read_snapshot(string path)
{
    double Db1,Db2,Db3;
    mvector Mv1;

    FILE * f1;

    f1 = fopen(path.c_str(),"r");

    //read position and rotation of rigid objects
    //M1
    fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
    KIN_M1_RB.pos.x = Db1;
    KIN_M1_RB.pos.y = Db2;
    KIN_M1_RB.pos.z = Db3;

    fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
    KIN_M1_RB.rot.x = Db1;
    KIN_M1_RB.rot.y = Db2;
    KIN_M1_RB.rot.z = Db3;

    //M2
    fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
    KIN_M2_RB.pos.x = Db1;
    KIN_M2_RB.pos.y = Db2;
    KIN_M2_RB.pos.z = Db3;

    fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
    KIN_M2_RB.rot.x = Db1;
    KIN_M2_RB.rot.y = Db2;
    KIN_M2_RB.rot.z = Db3;

    //CC-NL-Junction
    fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
    CCNLJ_RB.pos.x = Db1;
    CCNLJ_RB.pos.y = Db2;
    CCNLJ_RB.pos.z = Db3;

    fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
    CCNLJ_RB.rot.x = Db1;
    CCNLJ_RB.rot.y = Db2;
    CCNLJ_RB.rot.z = Db3;

    //cargo
    fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
    CAR_RB.pos.x = Db1;
    CAR_RB.pos.y = Db2;
    CAR_RB.pos.z = Db3;

    fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
    CAR_RB.rot.x = Db1;
    CAR_RB.rot.y = Db2;
    CAR_RB.rot.z = Db3;

    //read bead positions
    //M1
    for(int i=0;i<KIN_N_beads;i++)
    {
        fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
        KIN_M1_beads[i].x = Db1;
        KIN_M1_beads[i].y = Db2;
        KIN_M1_beads[i].z = Db3;
    }

    //M2
    for(int i=0;i<KIN_N_beads;i++)
    {
        fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
        KIN_M2_beads[i].x = Db1;
        KIN_M2_beads[i].y = Db2;
        KIN_M2_beads[i].z = Db3;
    }

    //CC
    for(int i=0;i<CC_N_beads;i++)
    {
        fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
        CC_beads[i].x = Db1;
        CC_beads[i].y = Db2;
        CC_beads[i].z = Db3;
    }

    //CC-NL-Junction
    for(int i=0;i<4;i++)
    {
        fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
        CCNLJ_beads[i].x = Db1;
        CCNLJ_beads[i].y = Db2;
        CCNLJ_beads[i].z = Db3;
    }

    //cargo
    for(int i=0;i<2;i++)
    {
        fscanf(f1,"%lf %lf %lf\n",&Db1,&Db2,&Db3);
        CAR_beads[i].x = Db1;
        CAR_beads[i].y = Db2;
        CAR_beads[i].z = Db3;
    }

    //angular data
    fscanf(f1,"%lf %lf\n",&Db1,&Db2);

    CAR_obs_ang = Db1;
    CAR_obs_ang_shift = Db2;

    fclose(f1);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
/////////////////////////////MISC/////////////////////////////

//defining a function to invert M1 and M2 roles
void KIN_invert_roles()
{
    //invert parameters
    eps_nl_M1 = eps_nl;
    eps_nl_M2 = 0.0;

    eps_tbs_M1 = eps_tbs_strong;
    eps_tbs_M2 = eps_tbs_weak;
}

//defining a function to create a vmd topology representation file for a cg molecule
void create_VMD_topo_KIN()
{
    FILE * f1;

    //topology file 1
    f1 = fopen(path_topo.c_str(),"w");

    fprintf(f1,"mol new %s waitfor all\n",path_str_out_M2.c_str());
    fprintf(f1,"mol new %s waitfor all\n",path_str_out_CC.c_str());
    fprintf(f1,"mol new %s waitfor all\n",path_str_out_CCNLJ.c_str());
    fprintf(f1,"mol new KIN-MT.xyz waitfor all\n");
    fprintf(f1,"mol new %s waitfor all\n\n",path_str_out_CAR.c_str());

    fprintf(f1,"source topo-M1\n");
    fprintf(f1,"source topo-M2\n");
    fprintf(f1,"source topo-CC\n");
    fprintf(f1,"source topo-CCNLJ\n");
    fprintf(f1,"source topo-MT\n");
    fprintf(f1,"source topo-CAR\n\n");

    fprintf(f1,"axes location off\n");
    fprintf(f1,"axes location off\n\n");

    fprintf(f1,"mol top 0\n");
    fprintf(f1,"animate goto start\n\n");

    fprintf(f1,"display culling on\n");
    fprintf(f1,"display depthcue off\n\n");

    fprintf(f1,"mol top 1\n");
    fprintf(f1,"display resetview\n");
    fprintf(f1,"mol top 0\n\n");

    fclose(f1);

    //topology file 2
    f1 = fopen("topo-M1","w");

    fprintf(f1,"mol top %d\n\n",0);

    fprintf(f1,"mol modstyle 0 %d VDW 1.0\n",0);
    fprintf(f1,"mol modcolor 0 %d ColorID 0\n",0);
    fprintf(f1,"mol modselect 0 %d index 0 to 6\n",0);
    fprintf(f1,"mol smoothrep %d 0 5\n\n",0);

    fprintf(f1,"mol addrep %d\n",0);
    fprintf(f1,"mol modstyle 1 %d VDW 1.0\n",0);
    fprintf(f1,"mol modcolor 1 %d ColorID 6\n",0);
    fprintf(f1,"mol modselect 1 %d index 7 to 25\n",0);
    fprintf(f1,"mol smoothrep %d 1 5\n\n",0);

    fprintf(f1,"mol addrep %d\n",0);
    fprintf(f1,"mol modstyle 2 %d VDW 1.0\n",0);
    fprintf(f1,"mol modcolor 2 %d ColorID 4\n",0);
    fprintf(f1,"mol modselect 2 %d index 26 to 39\n",0);
    fprintf(f1,"mol smoothrep %d 2 5\n\n",0);

    for(int i=0;i<KIN_N_beads;i++)
    {
        fprintf(f1,"set sel [atomselect top \"index %d\"]\n",i);
        fprintf(f1,"$sel set radius {%lf}\n",KIN_radii[i]);
    }

    fclose(f1);

    //topology file 3
    f1 = fopen("topo-M2","w");

    fprintf(f1,"mol top %d\n\n",1);

    fprintf(f1,"mol modstyle 0 %d VDW 1.0\n",1);
    fprintf(f1,"mol modcolor 0 %d ColorID 1\n",1);
    fprintf(f1,"mol modselect 0 %d index 0 to 6\n",1);
    fprintf(f1,"mol smoothrep %d 0 5\n\n",1);

    fprintf(f1,"mol addrep %d\n",1);
    fprintf(f1,"mol modstyle 1 %d VDW 1.0\n",1);
    fprintf(f1,"mol modcolor 1 %d ColorID 6\n",1);
    fprintf(f1,"mol modselect 1 %d index 7 to 25\n",1);
    fprintf(f1,"mol smoothrep %d 1 5\n\n",1);

    fprintf(f1,"mol addrep %d\n",1);
    fprintf(f1,"mol modstyle 2 %d VDW 1.0\n",1);
    fprintf(f1,"mol modcolor 2 %d ColorID 7\n",1);
    fprintf(f1,"mol modselect 2 %d index 26 to 39\n",1);
    fprintf(f1,"mol smoothrep %d 2 5\n\n",1);

    for(int i=0;i<KIN_N_beads;i++)
    {
        fprintf(f1,"set sel [atomselect top \"index %d\"]\n",i);
        fprintf(f1,"$sel set radius {%lf}\n",KIN_radii[i]);
    }

    fclose(f1);

    //topology file 4
    f1 = fopen("topo-CC","w");

    fprintf(f1,"mol top %d\n\n",2);

    fprintf(f1,"mol modstyle 0 %d VDW 1.0\n",2);
    fprintf(f1,"mol modcolor 0 %d ColorID 6\n",2);
    fprintf(f1,"mol modselect 0 %d index 0 to %d\n",2,CC_N_beads);
    fprintf(f1,"mol smoothrep %d 0 5\n\n",2);

    for(int i=0;i<CC_N_beads;i++)
    {
        fprintf(f1,"set sel [atomselect top \"index %d\"]\n",i);
        fprintf(f1,"$sel set radius {%lf}\n",CC_radii);
    }

    fclose(f1);

    //topology file 5
    f1 = fopen("topo-CCNLJ","w");

    fprintf(f1,"mol top %d\n\n",3);

    fprintf(f1,"mol modstyle 0 %d VDW 1.0\n",3);
    fprintf(f1,"mol modcolor 0 %d ColorID 6\n",3);
    fprintf(f1,"mol modselect 0 %d index 0\n",3);
    fprintf(f1,"mol smoothrep %d 0 5\n\n",3);

    fprintf(f1,"mol addrep %d\n",3);
    fprintf(f1,"mol modstyle 1 %d VDW 1.0\n",3);
    fprintf(f1,"mol modcolor 1 %d ColorID 4\n",3);
    fprintf(f1,"mol modselect 1 %d index 1\n",3);
    fprintf(f1,"mol smoothrep %d 1 5\n\n",3);

    fprintf(f1,"mol addrep %d\n",3);
    fprintf(f1,"mol modstyle 2 %d VDW 1.0\n",3);
    fprintf(f1,"mol modcolor 2 %d ColorID 7\n",3);
    fprintf(f1,"mol modselect 2 %d index 2\n",3);
    fprintf(f1,"mol smoothrep %d 2 5\n\n",3);

    fprintf(f1,"set sel [atomselect top \"index %d\"]\n",0);
    fprintf(f1,"$sel set radius {%lf}\n",CC_radii);

    fprintf(f1,"set sel [atomselect top \"index %d\"]\n",1);
    fprintf(f1,"$sel set radius {%lf}\n",1.9);

    fprintf(f1,"set sel [atomselect top \"index %d\"]\n",2);
    fprintf(f1,"$sel set radius {%lf}\n",1.9);

    fclose(f1);

    //topology file 6
    f1 = fopen("topo-MT","w");

    fprintf(f1,"mol top %d\n\n",4);

    fprintf(f1,"mol modstyle 0 %d VDW 1.5\n",4);
    fprintf(f1,"mol modcolor 0 %d ColorID 10\n",4);
    fprintf(f1,"mol modselect 0 %d index 0 to %d\n\n",4,MT_Nrep*13-1);

    fprintf(f1,"mol addrep %d\n",4);
    fprintf(f1,"mol modstyle 1 %d VDW 1.5\n",4);
    fprintf(f1,"mol modcolor 1 %d ColorID 23\n",4);
    fprintf(f1,"mol modselect 1 %d index %d to %d\n\n",4,MT_Nrep*13,2*MT_Nrep*13);

    for(int i=0;i<2*MT_Nrep*13;i++)
    {
        fprintf(f1,"set sel [atomselect top \"index %d\"]\n",i);
        fprintf(f1,"$sel set radius {%lf}\n",MT_proto_rad);
    }

    fclose(f1);

    //topology file 7
    f1 = fopen("topo-CAR","w");

    fprintf(f1,"mol top %d\n\n",5);

    fprintf(f1,"mol modstyle 0 %d VDW 1.0\n",5);
    fprintf(f1,"mol modcolor 0 %d ColorID 7\n",5);
    fprintf(f1,"mol modselect 0 %d index 0 to 1\n",5);
    fprintf(f1,"mol smoothrep %d 0 5\n\n",5);

    fprintf(f1,"set sel [atomselect top \"index 0\"]\n");
    fprintf(f1,"$sel set radius {%lf}\n",CAR_rad);

    fprintf(f1,"set sel [atomselect top \"index 1\"]\n");
    fprintf(f1,"$sel set radius {%lf}\n",CC_radii);

    fclose(f1);
}

/////////////////////////////END//////////////////////////////
//////////////////////////////////////////////////////////////

#endif
