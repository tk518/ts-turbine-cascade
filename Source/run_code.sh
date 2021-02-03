#!/bin/bash
#Make executable with chmod command

#Create the steady input
python make_design.py

# Run the steady code through Turbostream
mpirun -npernode 1 -np 1 turbostream input_1_psi_1.60_phi_0.50_Ma_0.70_slip.hdf5 output_1_psi_1.60_phi_0.50_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_1_psi_1.60_phi_0.70_Ma_0.70_slip.hdf5 output_1_psi_1.60_phi_0.70_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_1_psi_1.60_phi_0.90_Ma_0.70_slip.hdf5 output_1_psi_1.60_phi_0.90_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_1_psi_1.60_phi_1.10_Ma_0.70_slip.hdf5 output_1_psi_1.60_phi_1.10_Ma_0.70_slip 1

# Convert steady output to unsteady input
python convert_unsteady.py

# Run the unsteady code through Turbostream
mpirun -npernode 1 -np 1 turbostream input_2_psi_1.60_phi_0.50_Ma_0.70_slip.hdf5 output_2_psi_1.60_phi_0.50_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_2_psi_1.60_phi_0.70_Ma_0.70_slip.hdf5 output_2_psi_1.60_phi_0.70_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_2_psi_1.60_phi_0.90_Ma_0.70_slip.hdf5 output_2_psi_1.60_phi_0.90_Ma_0.70_slip 1

mpirun -npernode 1 -np 1 turbostream input_2_psi_1.60_phi_1.10_Ma_0.70_slip.hdf5 output_2_psi_1.60_phi_1.10_Ma_0.70_slip 1

# Delete the unwanted probe files
find -regex '.*_probe_[0-9]+.hdf5$' -delete
