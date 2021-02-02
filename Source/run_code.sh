#!/bin/bash
#Create the steady input
python make_design.py

'''
# Run the steady code through Turbostream
mpirun -npernode 1 -np 1 turbostream input_1_psi_1.60_phi_0.45_Ma_0.70_slip.hdf5 output_1_psi_1.60_phi_0.45_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_1_psi_1.60_phi_0.60_Ma_0.70_slip.hdf5 output_1_psi_1.60_phi_0.60_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_1_psi_1.60_phi_0.80_Ma_0.70_slip.hdf5 output_1_psi_1.60_phi_0.80_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_1_psi_1.60_phi_1.00_Ma_0.70_slip.hdf5 output_1_psi_1.60_phi_1.00_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_1_psi_1.60_phi_1.15_Ma_0.70_slip.hdf5 output_1_psi_1.60_phi_1.15_Ma_0.70_slip 1

# Convert steady output to unsteady input
$ python convert_unsteady.py

# Run the unsteady code through Turbostream
mpirun -npernode 1 -np 1 turbostream input_2_psi_1.60_phi_0.45_Ma_0.70_slip.hdf5 output_2_psi_1.60_phi_0.45_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_2_psi_1.60_phi_0.60_Ma_0.70_slip.hdf5 output_2_psi_1.60_phi_0.60_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_2_psi_1.60_phi_0.80_Ma_0.70_slip.hdf5 output_2_psi_1.60_phi_0.80_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_2_psi_1.60_phi_1.00_Ma_0.70_slip.hdf5 output_2_psi_1.60_phi_1.00_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_2_psi_1.60_phi_1.15_Ma_0.70_slip.hdf5 output_2_psi_1.60_phi_1.15_Ma_0.70_slip 1
'''