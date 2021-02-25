"""This script reads and plots the blade-to-blade probes."""
import numpy as np  # Multidimensional array library
import probe  # Code for reading TS probe output
import matplotlib.pyplot as plt  # Plotting library
from ts import ts_tstream_reader, ts_tstream_cut  # TS grid reader
from ts import ts_tstream_reader, ts_tstream_patch_kind
import os

#
# Set variables here
#

Phi = [0.4, 0.8, 1.2]
Psi = [1.6]
Ma = [0.7]
slip = False

for Psii in Psi:

    for Phii in Phi:   
        n = 0

        for Mai in Ma:
            if slip == True:
                output_file_name = 'output_2_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f_slip' % Mai  # Location of TS output file
            else:
                 output_file_name = 'output_2_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai  # Location of TS output file
            # Load the grid 
            tsr = ts_tstream_reader.TstreamReader()
            g = tsr.read(output_file_name + '.hdf5')

            
            #pid_probe = 8   # Patch ID of probes
            print('***')
            bid_probe = []
            for bid in g.get_block_ids():
                for pid in g.get_patch_ids(bid):
                    patch = g.get_patch(bid,pid)
                    if patch.kind == ts_tstream_patch_kind.probe:
                        rpm = g.get_bv('rpm',bid)
                        row_str = 'STATOR' if rpm==0 else 'ROTOR'
                        di = patch.ien- patch.ist 
                        dj = patch.jen- patch.jst
                        dk = patch.ken- patch.kst
                        if dj == 1:
                            bid_probe.append((bid,pid))
                        print('%s, bid=%d, pid=%d, di=%d, dj=%d, dk=%d' % (row_str, bid, pid, di, dj, dk))
            print('***')

            print 'bid_probe: ', bid_probe
            print 'length bid_probe: ', len(bid_probe)
            # This next section contains code to read in the data and process it into a
            # convenient form. Only a vague undestanding of this section is needed.
            # Store all probe patches in a list
            Dat = []
            for i in range(len(bid_probe)):
                bpi = bid_probe[i][0]
                ppi = bid_probe[i][1]

                print('Reading probe %d of %d' % (i+1, len(bid_probe)))
                print 'bpi: ', bpi
                print 'ppi: ', ppi
                # Determine the number of grid points on probe patches
                # (We index the TS grid using i = streamwise, j = spanwise, k =
                # pitchwise)
                p = g.get_patch(bpi,ppi)
                di = p.ien - p.ist
                dj = p.jen - p.jst
                dk = p.ken - p.kst
                probe_shape = [di, dj, dk]  # Numbers of points in i, j, k directions

                # Assemble file names for the probes using % substitution
                probe_name = (output_file_name + '_probe_%d_%d.dat' % (bpi,ppi))

                # Read the probes: The dictionary is keyed by variable name; the values
                # are numpy arrays with indexes 
                # [i = streamwise, j = spanwise, k = pitchwise, n = timewise] 
                Dat.append(probe.read_dat(probe_name, probe_shape))

            # Here we extract some parameters from the TS grid to use later
            rpm = g.get_bv('rpm',g.get_nb()-1)  # RPM in rotor row 
            print(rpm)
            Omega = rpm / 60. * np.pi * 2.
            print(Omega)
            cp = g.get_av('cp')  # Specific heat capacity at const p
            ga = g.get_av('ga')  # Specific heat ratio
            rgas = cp * (1.-1./ga)

            # Get information about time discretisation from TS grid
            freq = g.get_av('frequency')  # Blade passing frequency
            print(freq)
            ncycle = g.get_av('ncycle')  # Number of cycles
            nstep_cycle = g.get_av('nstep_cycle')  # Time steps per cycle
            print(nstep_cycle)
            nstep_save_probe = g.get_av('nstep_save_probe')  # Time steps per cycle
            # Individual time step in seconds = blade passing period / steps per cycle
            dt = 1./freq/float(nstep_cycle)

            # Number of time steps = num cycles * steps per cycle
            # nt = ncycle * nstep_cycle
            nt = np.shape(Dat[0]['ro'])[-1]
            print(dt)
            # Make non-dimensional time vector = time in seconds * blade passing frequency


            # Arbitrary datum temperatures for entropy level
            Pdat=16e5
            Tdat=1600.
            bid_rotor = g.get_nb()-1
            b_rotor = g.get_block(bid_rotor)

            # Get cuts at rotor inlet and exit to form Cp
            rotor_inlet = ts_tstream_cut.TstreamStructuredCut()
            rotor_inlet.read_from_grid(
                    g, Pdat, Tdat, bid_rotor,
                    ist = 0, ien=1,  # First streamwise
                    jst = 0, jen=b_rotor.nj,  # All radial
                    kst = 0, ken=b_rotor.nk  # All pitchwise
                    )
            rotor_outlet = ts_tstream_cut.TstreamStructuredCut()
            rotor_outlet.read_from_grid(
                    g, Pdat, Tdat, bid_rotor,
                    ist = b_rotor.ni-1, ien=b_rotor.ni,  # Last streamwise
                    jst = 0, jen=b_rotor.nj,  # All radial
                    kst = 0, ken=b_rotor.nk  # All pitchwise
                    )

            # Pressure references
            _, Po1 = rotor_inlet.mass_avg_1d('pstag')
            _, P1 = rotor_inlet.area_avg_1d('pstat')
            _, P2 = rotor_outlet.area_avg_1d('pstat')

            # Temperature references
            _, T1 = rotor_inlet.area_avg_1d('tstat')
            _, T2 = rotor_outlet.area_avg_1d('tstat')

            rpms = [g.get_bv('rpm',b) for b in g.get_block_ids()]
            print(rpms)
            # Get secondary vars, things like static pressure, rotor-relative Mach, etc.
            Dat = [probe.secondary(Di, rpmi, cp, ga, P1, T1) for Di, rpmi in zip(Dat,rpms)]

            # Determine sector size
            dtheta_sector = 2. * np.pi * g.get_bv('fracann',0)

            # Angular movement per time step
            dtheta_step = Omega * dt
            # We must add this as an offset between the blade rows because the mesh is
            # moved to the next time step before the probes are written out!

            # Finished reading data, now make some plots

            #Creating folder to store images
            path = os.getcwd()
            print(path)
            if slip == True:
                newpath = os.path.join(path, 'psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f_slip' % Mai)
            else:
                newpath = os.path.join(path, 'psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai)
            try:
                os.mkdir(newpath)
            except OSError:
                print ("Creation of the directory %s failed" % newpath)
            else:
                print ("Successfully created the directory %s " % newpath)

            '''
            for stepsize in range(nt): #(1, 97) gets all the steps
                f,a = plt.subplots()  # Create a figure and axis to plot into
                #plt.set_cmap('cubehelix_r')
                lev = np.linspace(-1.2,0.,61)
                probe.render_frame(a, Dat,'Cp', stepsize, lev, Omega, dt*nstep_save_probe,nstep_cycle/nstep_save_probe, dtheta_sector, Po1, P2)
                a.axis('equal')
                plt.grid(False)
                a.axis('off')
                a.set_ylim([0.,dtheta_sector*Dat[0]['r'][0,0,0,0]])
                plt.tight_layout()  # Remove extraneous white space
                plt.savefig('pressure_coeff_%d.png'%stepsize,dpi=200)
                plt.close(f)
            
            for stepsize in range(nt): #(1, 97) gets all the steps
                f,a = plt.subplots()  # Create a figure and axis to plot into
                #plt.set_cmap('cubehelix_r')
                lev = np.linspace(-0.015,0.015,201)
                probe.render_frame(a, Dat,'pmean', stepsize, lev, Omega, dt*nstep_save_probe,nstep_cycle/nstep_save_probe, dtheta_sector, Po1, P2)
                a.axis('equal')
                plt.grid(False)
                a.axis('off')
                a.set_ylim([0.,dtheta_sector*Dat[0]['r'][0,0,0,0]])
                plt.tight_layout()  # Remove extraneous white space
                if slip == True:
                    plt.savefig(newpath + '/mean_static_pressure_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f_slip' % Mai +'_%d.png' %stepsize, dpi=200)
                else:
                    plt.savefig(newpath + '/mean_static_pressure_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai +'_%d.png' %stepsize, dpi=200)
                plt.close(f)


            # Entropy
            for stepsize in range(nt): #(1, 97) gets all the steps
                f,a = plt.subplots()  # Create a figure and axis to plot into
                #plt.set_cmap('cubehelix_r')
                lev = np.linspace(-8.,35.0,21)
                probe.render_frame(a, Dat,'ds', stepsize, lev, Omega, dt*nstep_save_probe,nstep_cycle/nstep_save_probe, dtheta_sector, Po1, P2)
                a.axis('equal')
                plt.grid(False)
                a.axis('off')
                a.set_ylim([0.,dtheta_sector*Dat[0]['r'][0,0,0,0]])
                plt.tight_layout()  # Remove extraneous white space
                if slip == True:
                    plt.savefig(newpath + '/Entropy_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f_slip' % Mai + '_%d.png' %stepsize, dpi=200)
                else:
                    plt.savefig(newpath + '/Entropy_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai + '_%d.png' %stepsize, dpi=200)
                plt.close(f)
            '''

            #Entropy contour with pressure plot
            for stepsize in range(nt): #(1, 97) gets all the steps
                f,a = plt.subplots()  # Create a figure and axis to plot into
                #plt.set_cmap('cubehelix_r')
                lev = np.linspace(-0.015,0.015,201)
                probe.render_frame(a, Dat,'pmean', stepsize, lev, Omega, dt*nstep_save_probe,nstep_cycle/nstep_save_probe, dtheta_sector, Po1, P2)
                lev = np.linspace(-8.,35.0,21)
                probe.render_frame(a, Dat,'ds', stepsize, lev, Omega, dt*nstep_save_probe,nstep_cycle/nstep_save_probe, dtheta_sector, Po1, P2)
                a.axis('equal')
                plt.grid(False)
                a.axis('off')
                a.set_ylim([0.,dtheta_sector*Dat[0]['r'][0,0,0,0]])
                plt.tight_layout()  # Remove extraneous white space
                if slip == True:
                    plt.savefig(newpath + '/combined_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f_slip' % Mai +'_%d.png' %stepsize, dpi=200)
                else:
                    plt.savefig(newpath + '/combined_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai +'_%d.png' %stepsize, dpi=200)
                plt.close(f)

            #plt.show()  # Render the plot