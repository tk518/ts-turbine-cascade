#!/bin/bash
#Make executable with chmod command

#Create the steady input
python make_design.py

# Run the steady code through Turbostream
for i in 0.80 1.20 1.60 2.00 2.40
do 
    for x in 0.40 0.60 0.80 1.00 1.20
    do
        echo "mpirun -npernode 1 -np 1 turbostream input_1_psi_"$i"_phi_"$x"_Ma_0.70.hdf5 output_1_psi_"$i"_phi_"$x"_Ma_0.70 1"
        mpirun -npernode 1 -np 1 turbostream input_1_psi_"$i"_phi_"$x"_Ma_0.70_slip.hdf5 output_1_psi_"$i"_phi_"$x"_Ma_0.70_slip 1
    done
done
# Convert steady output to unsteady input
python convert_unsteady.py

# Run the unsteady code through Turbostream
for i in 0.80 1.20 1.60 2.00 2.40
do 
    for x in 0.40 0.60 0.80 1.00 1.20
    do
        mpirun -npernode 1 -np 1 turbostream input_2_psi_"$i"_phi_"$x"_Ma_0.70_slip.hdf5 output_2_psi_"$i"_phi_"$x"_Ma_0.70_slip 1
    done
done

# Delete the unwanted probe files
find -regex '.*_probe_[0-9]+.hdf5$' -delete
