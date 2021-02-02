"""
Choose parameters for a turbine stage and generate a corresponding CFD model.
"""
import design

fname = 'input_1'  # Desired name of the TS input hdf5
# phi = 0.8  # Flow coefficient (0.4 to 1.2)
# psi = 1.6  # Stage loading coefficient (0.8 to 2.4)
Lam = 0.5  # Degree of reaction (0.4 to 0.6)
# Ma = 0.9  # Vane exit Mach number (0.6 to 0.9)
eta = 0.9  # Polytropic efficiency (leave this for now)
gap_chord = 0.5  # Spacing between stator and rotor
slip_vane = True  # Apply non-slip condition to vane surface
guess_file = None #'guess.hdf5'  # Solution to use as initial guess, or None


# Check that the code works for many Mach
if slip_vane == False:
    for Psii in [1.6]:

        for Phii in [0.45, 0.6, 0.8, 1.0, 1.15]:

            for Mai in [0.7]:

                # Create a file name for this Mach using % substitution
                fname_now = fname + '_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai + '.hdf5'

                # Call out to the design generation code
                # It is complicated so best to think of it as a black box!
                design.generate(fname_now, phi, psi, Lam, Ma, eta, gap_chord, slip_vane, guess_file )
    Else:
        for Psii in [1.6]:

            for Phii in [0.45, 0.6, 0.8, 1.0, 1.15]:

                for Mai in [0.7]:

                    # Create a file name for this Mach using % substitution
                    fname_now = fname + '_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai + '_slip.hdf5'

                    # Call out to the design generation code
                    # It is complicated so best to think of it as a black box!
                    design.generate(fname_now, phi, psi, Lam, Ma, eta, gap_chord, slip_vane, guess_file )