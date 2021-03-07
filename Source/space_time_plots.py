"""This script reads in an unsteady solution and runs the simple hole model."""
import numpy as np  # Multidimensional array library
import probe  # Code for reading TS probe output
import model  # Simple hole model
import matplotlib.pyplot as plt  # Plotting library
from matplotlib.pyplot import cm
from ts import ts_tstream_reader  # TS grid reader
from ts import ts_tstream_cut  # TS cutter

#
# Set variables here
#

def rms(x):
        ms = np.sqrt(np.mean(x**2, axis = 0))
        return(ms)

Data={}

Phi = [0.40]
Psi = [1.6]
Ma =  [0.70]

slip = False

if slip == False:
        for Psii in Psi:

                for Phii in Phi:
                        n = 0
                        #Mach = []

                        for Mai in Ma:

                                output_file_name = 'output_2_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai  # Location of TS output file

                                # We identify a region of the grid using block and patch IDs
                                pid_probe_ps = 9  # Patch ID of surface probe on pressure side
                                pid_probe_ss = 10  # Patch ID of surface probe on suction side
                                pid_probe_ps_free = 11  # Patch ID of free-stream probe on pressure side
                                pid_probe_ss_free = 12  # Patch ID of free-stream probe on suction side

                                #
                                # This next section contains code to read in the data and process it into a
                                # convenient form. Only a vague undestanding of this section is needed.
                                #

                                # Load the grid 
                                tsr = ts_tstream_reader.TstreamReader()
                                g = tsr.read(output_file_name + '.hdf5')

                                # Determine number of blades in each row
                                bids = [0,g.get_nb()-1]
                                fracann = np.array([g.get_bv('fracann',bi) for bi in bids])
                                nblade = np.array([g.get_bv('nblade',bi) for bi in bids])
                                nb_row = np.round(fracann * nblade)
                                bid_probe = int(nb_row[0])  # Block ID where probes are located


                                # Determine the number of grid points on probe patches
                                # (We index the TS grid using i = streamwise, j = spanwise, k = pitchwise)
                                p = g.get_patch(bid_probe,pid_probe_ps)
                                di = p.ien - p.ist
                                dj = p.jen - p.jst
                                probe_shape = [di, dj, 1]  # Numbers of points in i, j, k directions

                                # Index for the mid-span
                                jmid = int(dj/2)

                                # Assemble file names for the probes using % substitution
                                probe_name_ps = output_file_name + '_probe_%d_%d.dat' % (bid_probe,pid_probe_ps)
                                probe_name_ss = output_file_name + '_probe_%d_%d.dat' % (bid_probe,pid_probe_ss)
                                probe_name_ps_free = output_file_name + '_probe_%d_%d.dat' % (bid_probe,pid_probe_ps_free)
                                probe_name_ss_free = output_file_name + '_probe_%d_%d.dat' % (bid_probe,pid_probe_ss_free)

                                # Read the probes
                                # The probe data are separate dictionary for each surface of the blade The
                                # dictionary is keyed by variable name; the values are numpy arrays with
                                # indexes [i = streamwise, j = spanwise, k = pitchwise, n = timewise]
                                # For example, to get density at time instant n as a function of
                                # axial distance at mid-radius:
                                #   Dat_ps['ro'][:,jmid,0,n]
                                Dat_ps = probe.read_dat(probe_name_ps, probe_shape)
                                Dat_ss = probe.read_dat(probe_name_ss, probe_shape)
                                Dat_ps_free = probe.read_dat(probe_name_ps_free, probe_shape)
                                Dat_ss_free = probe.read_dat(probe_name_ss_free, probe_shape)

                                # Here we extract some parameters from the TS grid to use later
                                rpm = g.get_bv('rpm',1)  # RPM in rotor row
                                cp = g.get_av('cp')  # Specific heat capacity at const p
                                ga = g.get_av('ga')  # Specific heat ratio

                                # Get information about time discretisation from TS grid
                                freq = g.get_av('frequency')  # Blade passing frequency
                                ncycle = g.get_av('ncycle')  # Number of cycles
                                nstep_cycle = g.get_av('nstep_cycle')  # Time steps per cycle
                                nstep_save_probe = g.get_av('nstep_save_probe')  # Time steps per cycle
                                # Individual time step in seconds = blade passing period / steps per cycle
                                dt = 1./freq/float(nstep_cycle)*float(nstep_save_probe)
                                # Number of time steps = num cycles * steps per cycle
                                # nt = ncycle * nstep_cycle
                                nt = np.shape(Dat_ps['ro'])[-1]
                                print(nt)
                                # Make non-dimensional time vector = time in seconds * blade passing frequency
                                ft = np.linspace(0.,float(nt-1)*dt,nt) * freq

                                # Get secondary vars, things like static pressure, rotor-relative Mach, etc.
                                Dat_ps = probe.secondary(Dat_ps, rpm, cp, ga, 1, 1)
                                Dat_ss = probe.secondary(Dat_ss, rpm, cp, ga, 1, 1)
                                Dat_ps_free = probe.secondary(Dat_ps_free, rpm, cp, ga, 1, 1)
                                Dat_ss_free = probe.secondary(Dat_ss_free, rpm, cp, ga, 1, 1)

                                # Cut the rotor inlet
                                Pdat = 1e5
                                Tdat = 300.
                                b = g.get_block(bid_probe)
                                rotor_inlet = ts_tstream_cut.TstreamStructuredCut()
                                rotor_inlet.read_from_grid(
                                        g, Pdat, Tdat, bid_probe,
                                        ist = 0, ien=1,  # First streamwise
                                        jst = 0, jen=b.nj,  # All radial
                                        kst = 0, ken=b.nk # All pitchwise
                                        )

                                # Get mass averaged rotor inlet relative stagnation conditions
                                _, Po1 = rotor_inlet.mass_avg_1d('pstag_rel')
                                _, To1 = rotor_inlet.mass_avg_1d('tstag_rel')

                                #
                                # Set up the simple hole model
                                #

                                # Choose a constant "pressure margin", or percentage increase of coolant
                                # relative to the inlet stagnation condition
                                PM = 0.01
                                Poc = (1. + PM) * Po1

                                # Fix a stagnation temperature ratio, i.e. the coolant is this much colder than
                                # the main-stream inlet stagnation condition
                                TR = 0.5
                                Toc = TR * To1

                                x = Dat_ps['x'][:,jmid,0,0]
                                #print(x.min())
                                #print(x.max())

                                # Pull out data for model

                                Pps = Dat_ps_free['pstat'][:,jmid,0,:]
                                Pss = Dat_ss_free['pstat'][:,jmid,0,:]
                                #
                                # Finished reading data, now make some plots

                                #key in form 'Ma_0.70_psi_1.60_phi_0.45'
                                Data['Ma_'+"{:.2f}".format(Mai)+'_psi_'+"{:.2f}".format(Psii)+'_phi_'+"{:.2f}".format(Phii)] = [Pps, Pss]
                                #find out which way round

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

                                f,a = plt.subplots()  # Create a figure and axis to plot into
                                P_stat_ps = Data['Ma_'+"{:.2f}".format(Mai)+'_psi_'+"{:.2f}".format(Psii)+'_phi_'+"{:.2f}".format(Phii)][0]
                                P_stat_ss = Data['Ma_'+"{:.2f}".format(Mai)+'_psi_'+"{:.2f}".format(Psii)+'_phi_'+"{:.2f}".format(Phii)][1]
                                lev = np.linspace(-0.015,0.015,21)
                                for i in range(len(x)):                     #Range of x values
                                    for ii in range(len(P_stat_ps[1]))      #Rage of t values
                                        a.contourf(-x[i], ii, P_stat_ps[i][ii], lev)
                                        a.contourf(x[i], ii, P_stat_ss[i][ii], lev)

                                a.axis('equal')
                                plt.grid(False)
                                a.axis('off')
                                plt.tight_layout()
                                plt.savefig(newpath + '/space_time_pressure_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai +'_%d.png' %stepsize, dpi=200)