"""
Choose parameters for a turbine stage and generate a corresponding CFD model.
"""
import design

for x in range(1, 10):   
    phi = 0.8  # Flow coefficient (0.4 to 1.2)
    psi = 1.6  # Stage loading coefficient (0.8 to 2.4)
    Lam = 0.5  # Degree of reaction (0.4 to 0.6)
    Ma = 0.1*x  # Vane exit Mach number (0.6 to 0.9)
    eta = 0.9  # Polytropic efficiency (leave this for now)
    gap_chord = 0.5  # Spacing between stator and rotor
    fname = 'input_Mach' +str(x)+ '.hdf5'  # Desired name of the TS input hdf5
    design.generate(fname, phi, psi, Lam, Ma, eta, gap_chord )
    
# Call out to the design generation code
# It is complicated so best to think of it as a black box!

