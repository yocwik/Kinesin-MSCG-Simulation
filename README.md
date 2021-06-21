# Readme File
This is a description of the code for simulating the motion of the kinesin motor using a multi-scale coarse grained model. The details of the model can be found in this [paper](https://www.biorxiv.org/content/10.1101/2021.04.05.438476v1.abstract).

## General Usage Notes
This code was written in C++ and can in theory be compiled and run on any system, given that all the mandatory files are presents and that an appropriate compiler is available. For compilation purpose, please compile the file main.cpp to produce the desired executable file.

## Usage Instructions
To run the simulation, simply execute the executable file generated by the compilation. As a default, this file will be referred to from here on as main.exe. This file can be provided with the following flags and inputs:
1. -s [int]: This flag, followed by an integer, provides a seed number for the random number generator in the simulation. It also adds the provided number to the names of all the output files of the simulation.
2. -x: Trajectory files will not be saved.
3. -f [float]: An external force whose magnitude in pN is the provide float argument will be applied to the cargo.
4. -fv [float]: An increasing external force will be applied to the cargo. The force increasing rate is the provided float argument in pN per time step.
5. -i: Inverse the roles of motor domains 1 and 2.
6. -d: Forced detachment of motor domain 1 from the microtubule. If -i was provided as an argument, motor domain 2 is the one to detach.
7. -a [float]: Accelerate the cargo dynamics by a factor of the provided float argument.
8. -r: Record a snap-shot of the system at the end of the simulation. The snap-shot is saved in a file called snapshot-conf.dat.
9. -l: Use the snap-shot file as the initial configuration of the system instead of init-conf.dat. Important note: if the -s flag was provided to generate the snap-shot, the same seed number has to be provided with the -s flag again in order for the simulation to work.
10. -e: Terminate the simulation when the stepping motor domain reaches a target binding site on the microtubule.

To change the simulation parameters, simply change the values of the numbers in the kin-param file. The details of the units of the different parameters can be found in the KIN_read_param function in the kinesin.h source code file. The frames value determines how many frames in total will be saved, which in turn determines for how long will the simulation run (unless the -e flag was provided). The save-period value determines how many time steps pass before a frame is saved.

## Mandatory Files
These are the files needed to compile and run the simulations.
Code files:
* main.cpp
* bd.h
* forces.h
* kinesin.h
* kinstr.h
* lists.h
* maths.h
* observables.h
* ranfun.h
* rbd.h
* rbdcg.h

Mandatory input files:
* CGmodel-Kin-coord.txt
* CGmodel-Kin-radii.txt
* DiffTen-2KIN.dat
* DiffTenChol-2KIN.dat
* init-conf.dat
* kin-param

Additional files:
* KIN-MT.xyz

## Output Files
By default the simulation will produce the following files as output. Note: if the -s flag is provided, some of these file names will have the seed number added to them.
Trajectory files:
* KIN-CAR.xyz
* KIN-CC.xyz
* KIN-CCNLJ.xyz
* KIN-M1.xyz
* KIN-M2.xyz

Observables file:
* obs-data.dat

Visualization scripts for VMD:
* topo
* topo-CAR
* topo-CC
* topo-CCNLJ
* topo-M1
* topo-M2
* topo-MT

## Observables
The details of the observables can be found in the observables.h source code file

## Visualization
The visualization scripts are meant to be used in VMD, which can be downloaded freely at [this link](http://www.ks.uiuc.edu/Research/vmd/).

In order to visualize the simulation trajectory, open the file KIN-M1.xyz in VMD and run the topo script using the source command.
