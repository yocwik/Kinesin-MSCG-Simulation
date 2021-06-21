#include <iostream>
#include <cstring>
#include <vector>
#include <ctime>

#include "maths.h"
#include "ranfun.h"
#include "bd.h"
#include "rbd.h"
#include "rbdcg.h"
#include "forces.h"
#include "lists.h"
#include "kinstr.h"
#include "kinesin.h"
#include "observables.h"

using namespace std;

int main(int argc,char * argv[])
{
    //should structure frames be saved
    bool write_str = true;

    //external force
    double ex_force = 0.0;

    //external force loading rate (pN per time step)
    double ex_force_inc = 0.0;

    //is external force variable
    bool ex_force_var = false;

    //should M1 and M2 roles be inverted
    bool inv_roles = false;

    //should TH detach from initial binding site and not rebind
    bool detach_TH = false;

    //should cargo dynamics be accelerated
    bool acc_cargo = false;

    //cargo dynamics acceleration factor
    double acc_cargo_factor = 1.0;

    //record configuration snapshot at the end of the simulation
    bool record_snapshot = false;

    //load snapshot as initial configuration
    bool load_snapshot = false;

    //end simulation when stepping motor reaches TBS
    bool terminate_at_TBS = false;

    //termination countdown variable
    int TBS_countdown = 10;

    //handling command line parameters
    for(int i=0;i<argc;i++)
    {
        if(argv[i][0] == '-')
        {
            switch(argv[i][1])
            {
                //seed number provided through command line argument
                case 's':
                    init += atoi(argv[i+1]);

                    //update structure paths
                    path_str_out_M1 = "KIN-M1-";
                    path_str_out_M2 = "KIN-M2-";
                    path_str_out_CC = "KIN-CC-";
                    path_str_out_CCNLJ = "KIN-CCNLJ-";
                    path_str_out_CAR = "KIN-CAR-";

                    path_str_out_M1.append(argv[i+1]);
                    path_str_out_M2.append(argv[i+1]);
                    path_str_out_CC.append(argv[i+1]);
                    path_str_out_CCNLJ.append(argv[i+1]);
                    path_str_out_CAR.append(argv[i+1]);

                    path_str_out_M1.append(".xyz");
                    path_str_out_M2.append(".xyz");
                    path_str_out_CC.append(".xyz");
                    path_str_out_CCNLJ.append(".xyz");
                    path_str_out_CAR.append(".xyz");

                    //update topology paths
                    path_topo = "topo-";
                    path_topo.append(argv[i+1]);

                    //update observables path
                    path_obs.append("-");
                    path_obs.append(argv[i+1]);

                    //update snapshot paths
                    //snapshot output paths
                    snapshot_path = "snapshot-conf-";

                    snapshot_path.append(argv[i+1]);

                    snapshot_path.append(".dat");

                    cout << "Seed set to: " << atoi(argv[i+1]) << endl;
                    break;

                //structure trajectories will no be saved
                case 'x':
                    write_str = false;

                    cout << "Structure trajectory files will not be saved." << endl;
                    break;

                //external force in pN is provided through command line argument
                case 'f':
                    //determine if force increases during simulation or is constant
                    if(argv[i][2]=='v')
                    {
                        ex_force_inc = double(atof(argv[i+1]));
                        ex_force_var = true;

                        cout << "External forces loading rate set to: " << ex_force_inc << " pN per time step" << endl;
                    } else {
                        ex_force = double(atof(argv[i+1]));

                        cout << "External forces set to: " << ex_force << " pN" << endl;
                    }

                    break;

                //inverse the roles of the M1 and M2
                case 'i':
                    inv_roles = true;

                    cout << "M1 and M2 roles are inverted." << endl;
                    break;

                //force TH detachment
                case 'd':
                    detach_TH = true;

                    cout << "TH is forced to detach from its initial binding site." << endl;
                    break;

                //accelerate cargo dynamics
                case 'a':
                    acc_cargo = true;
                    acc_cargo_factor = double(atof(argv[i+1]));

                    cout << "Cargo dynamics acceleration factor set to: " << acc_cargo_factor << endl;
                    cout << "Motor detachment is prevented." << endl;
                    break;

                //record snapshot at the end of simulation
                case 'r':
                    record_snapshot = true;

                    cout << "A snapshot of the system configuration will be saved at the end of the simulation." << endl;
                    break;

                //load snapshot as initial configuration
                case 'l':
                    load_snapshot = true;

                    cout << "Snapshot will be used as the initial configuration of the system." << endl;
                    break;

                //end simulation when TBS is reached
                case 'e':
                    terminate_at_TBS = true;

                    cout << "Simulation will terminate when the stepping motor reaches the TBS." << endl;
                    break;
            }
        }
    }

    path_obs.append(".dat");

    //modify initial configuration paths to snapshot paths
    if(load_snapshot)
    {
        init_conf_path = snapshot_path;
    }

    //initializing kinesin objects
    KIN_initialize_objects();

    //read simulation parameters from file
    KIN_read_param();

    //initializing coiled coil objects
    CC_initialize_objects();

    //creating coiled coil structure
    CC_create_CC_structure();

    //initializing rigid body dynamics parameters
    RBD_initialize_param();

    //initializing interaction lists
    KIN_initialize_lists();

    //initializing position and displacement lists
    KIN_create_pos_record();

    //initializing MT objects
    MT_initialize_objects(KIN_M1_beads);

    //creating MT xyz file
    //MT_create_xyz_file();
    //MT_create_xyz_file_mod_rad(-10.0);

    //initializing force parameters
    set_force_param();
    set_force_coeff();

    //initializing cargo objects
    CAR_initialize_objects();

    //creating cargo structure
    CAR_create_CAR_structure();

    //obtain initial value of observed angle
    CAR_obs_ang = CAR_get_obs_ang();

    //read initial system configuration from file
    read_snapshot(init_conf_path);

    //invert M1 and M2 roles
    if(inv_roles)
    {
        KIN_invert_roles();
    }

    //accelerate cargo dynamics
    if(acc_cargo)
    {
        //multiply cargo diffusion tensor by acceleration factor
        CAR_DT[0][0] *= acc_cargo_factor;
        CAR_DT[1][1] *= acc_cargo_factor;
        CAR_DT[2][2] *= acc_cargo_factor;
        CAR_DT[3][3] *= acc_cargo_factor;
        CAR_DT[4][4] *= acc_cargo_factor;
        CAR_DT[5][5] *= acc_cargo_factor;

        CAR_DT_Chol[0][0] *= sqrt(acc_cargo_factor);
        CAR_DT_Chol[1][1] *= sqrt(acc_cargo_factor);
        CAR_DT_Chol[2][2] *= sqrt(acc_cargo_factor);
        CAR_DT_Chol[3][3] *= sqrt(acc_cargo_factor);
        CAR_DT_Chol[4][4] *= sqrt(acc_cargo_factor);
        CAR_DT_Chol[5][5] *= sqrt(acc_cargo_factor);

        //make both motors tightly bound to MT
        eps_tbs_M1 = eps_tbs_strong;
		eps_tbs_M2 = eps_tbs_strong;

        detach_TH = false;
    }

    //initializing BAOAB
    KIN_initialize_BAOAB();

    //write initial frames
    if(write_str)
    {
        write_cg_beads_xyz(KIN_M1_beads,path_str_out_M1);
        write_cg_beads_xyz(KIN_M2_beads,path_str_out_M2);
        write_cg_beads_xyz(CC_beads,path_str_out_CC);
        write_cg_beads_xyz(CCNLJ_beads,path_str_out_CCNLJ);
        write_cg_beads_xyz(CAR_beads,path_str_out_CAR);
    }

    //create VMD topology file
    create_VMD_topo_KIN();

    //initialize observables file
    obs_write_file();

    //initial positions
    obs_init_M1_pos = KIN_M1_RB.pos;
    obs_init_M2_pos = KIN_M2_RB.pos;

    //reset frame and steps count
    sim_frame_count = 1;
    sim_step_count = 1;

    //start simulation
    while(sim_frame_count<=sim_N_frames)
    {
        //updating neighbor lists
        if(NLupdate)
        {
            KIN_update_neighbor_lists();
            KIN_reset_displacement_record();
            KIN_reset_pos_record();
        }

        //reseting forces
        KIN_M1_RB_force = set_vector(0,0,0);
        KIN_M2_RB_force = set_vector(0,0,0);
        CCNLJ_RB_force = set_vector(0,0,0);
        CAR_RB_force = set_vector(0,0,0);

        KIN_M1_RB_torque = set_vector(0,0,0);
        KIN_M2_RB_torque = set_vector(0,0,0);
        CCNLJ_RB_torque = set_vector(0,0,0);
        CAR_RB_torque = set_vector(0,0,0);

        reset_forces(KIN_M1_force);
        reset_forces(KIN_M2_force);
        reset_forces(KIN_M1_dock_cb_force);
        reset_forces(KIN_M1_dock_nl_force);
        reset_forces(KIN_M2_dock_cb_force);
        reset_forces(KIN_M2_dock_nl_force);
        reset_forces(CC_force);
        reset_forces(CCNLJ_force);
        reset_forces(CAR_force);

        //applying connectivity forces
        KIN_apply_forces_connectivity();

        //applying excluded volume forces
        KIN_apply_forces_excluded_volume_internal();
        KIN_apply_forces_excluded_volume_external();

        //apply angular forces
        CC_apply_forces_angular();

        //apply docking forces
        KIN_apply_forces_docking();

        //apply neck linker - CC junction forces
        KIN_apply_forces_CCNLJ();

        //update TBS indexes
        MT_TBS_M1_ind = MT_get_TBS_index(KIN_M1_RB.pos);
        MT_TBS_M2_ind = MT_get_TBS_index(KIN_M2_RB.pos);

        //update TBS positions
        MT_update_TBS_pos(MT_TBS_M1_pos,MT_TBS_M1_ind);
        MT_update_TBS_pos(MT_TBS_M2_pos,MT_TBS_M2_ind);

        //apply tbs forces
        if(!inv_roles)
        {
            //non stepping motor
            MT_apply_forces_TBS_M2();

            //stepping motor
            if(detach_TH)
            {
                if(MT_TBS_M1_ind.a==MT_TBS_M2_ind.a&&MT_TBS_M1_ind.b==MT_TBS_M2_ind.b+1)
                {
                    MT_apply_forces_TBS_M1();
                }
            } else {
                MT_apply_forces_TBS_M1();
            }
        } else {
            //non stepping motor
            MT_apply_forces_TBS_M1();

            //stepping motor
            if(detach_TH)
            {
                if(MT_TBS_M2_ind.a==MT_TBS_M1_ind.a&&MT_TBS_M2_ind.b==MT_TBS_M1_ind.b+1)
                {
                    MT_apply_forces_TBS_M2();
                }
            } else {
                MT_apply_forces_TBS_M2();
            }
        }

        //apply MT excluded volume forces
        MT_apply_forces_excluded_volume();

        //apply cargo forces
        CAR_apply_forces();

        //stall force
        CAR_force[0] += set_vector(ex_force,0,0);

        //resolve forces and torques of rigid bodies
        KIN_resolve_forces_torques();

        //get BAOAB fluctuations
        KIN_BAOAB_get_fluc();

        //integrate BAOAB time step
        KIN_BAOAB_integrate();

        //updating position and displacement records
        KIN_update_pos_record();
        KIN_update_displacement_record();

        //deciding if neighbor lists need updating
        NLupdate = update_neighbor_list(KIN_dr);

        //apply loading rate if applicable
        if(ex_force_var)
        {
            ex_force += ex_force_inc;
        }

        if(sim_step_count==sim_save_period)
        {
            cout << sim_frame_count << "/" << sim_N_frames << endl;

            if(write_str)
            {
                write_cg_beads_xyz_frame(KIN_M1_beads,path_str_out_M1);
                write_cg_beads_xyz_frame(KIN_M2_beads,path_str_out_M2);
                write_cg_beads_xyz_frame(CC_beads,path_str_out_CC);
                write_cg_beads_xyz_frame(CCNLJ_beads,path_str_out_CCNLJ);
                write_cg_beads_xyz_frame(CAR_beads,path_str_out_CAR);
            }

            obs_write_frame();

            //TBS termination conditions
            if(terminate_at_TBS)
            {
                if(!inv_roles)
                {
                    if((KIN_M1_RB.pos-MT_get_TBS_RB_pos(MT_TBS_M1_ind)).length()<=3.8)
                    {
                        TBS_countdown -= 1;
                    } else{
                        TBS_countdown = 10;
                    }
                } else {
                    if((KIN_M2_RB.pos-MT_get_TBS_RB_pos(MT_TBS_M2_ind)).length()<=3.8)
                    {
                        TBS_countdown -= 1;
                    } else{
                        TBS_countdown = 10;
                    }
                }

                if(TBS_countdown==0)
                {
                    break;
                }
            }

            //count frame
            sim_frame_count++;

            //step count reset
            sim_step_count = 1;

        } else {
            //count step
            sim_step_count++;
        }
    }

    //writing final configuration to snapshot file
    if(record_snapshot)
    {
        write_snapshot(snapshot_path);
    }

    return 0;
}
