"""This script reads in an unsteady solution and runs the simple hole model."""
import numpy as np  # Multidimensional array library
import probe  # Code for reading TS probe output
import model  # Simple hole model
import matplotlib.pyplot as plt  # Plotting library
from ts import ts_tstream_reader  # TS grid reader
from ts import ts_tstream_cut  # TS cutter

#
# Set variables here
#

def rms(x):
        ms = np.sqrt(np.mean(x**2))
        return(ms)
Data={}
Phi = [0.40, 0.60, 0.80, 1.00, 1.20]
Psi = [0.80, 1.20, 1.60, 2.00, 2.40]
Ma =  [0.70]

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
                g = tsr.read(output_file_name + '.hdf5', read_yplus=True)

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

                # Choose a hole position - 98 i positions on both sides
                ihole_ps = [20,40,60,80]
                ihole_ss = [20,40,60,80]

                # Pull out data for model
                roinf = np.stack((Dat_ps_free['ro'][ihole_ps,jmid,0,:],
                                Dat_ss_free['ro'][ihole_ss,jmid,0,:]))
                Vinf = np.stack((Dat_ps_free['vrel'][ihole_ps,jmid,0,:],
                                Dat_ss_free['vrel'][ihole_ss,jmid,0,:]))
                Pinf = np.stack((Dat_ps_free['pstat'][ihole_ps,jmid,0,:],
                                Dat_ss_free['pstat'][ihole_ss,jmid,0,:]))

                # Nondimensionalise data
                Pinf_Poc, roVinf_Po_cpToc = model.normalise(Poc, Toc, Pinf, roinf, Vinf, cp)

                # Assume constant Cd
                Cd = 0.7

                # Calculate BR
                BR = model.evaluate( Pinf_Poc, roVinf_Po_cpToc, Cd, ga )
                print 'BR: ', BR
                # Finished reading data, now make some plots

                # Plot the hole position
                if n == 0:
                        f,a = plt.subplots()  # Create a figure and axis to plot into
                        x = Dat_ps['x'][:,jmid,0,0]
                        rt_ps = Dat_ps['rt'][:,jmid,0,0]
                        rt_ss = Dat_ss['rt'][:,jmid,0,0]
                        a.plot(x,rt_ps,'-k')  # Blade pressure surface
                        a.plot(x,rt_ss,'-k')  # Blade suction surface
                        a.plot(x[ihole_ss],rt_ss[ihole_ss],'b*')  # SS hole location
                        a.plot(x[ihole_ps],rt_ps[ihole_ps],'r*')  # PS hole location 
                        plt.axis('equal')
                        plt.axis('off')
                        plt.tight_layout()  # Remove extraneous white space
                        plt.savefig('hole_posn_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '.pdf')  # Write out a pdf file

                #peak to peak amplitude

                #pressure side, if pressure side is [0]
                ptp_ps = np.max(BR.T[0],axis = 1)- np.min(BR.T[0],axis = 1)
                #suction side, if pressure side is [1]
                ptp_ss = np.max(BR.T[1],axis = 1)- np.min(BR.T[1],axis = 1)

                '''
                #rms BR
                #pressure side, if pressure side is [0]
                rms_ps = rms(BR.T[0])
                #suction side, if pressure side is [1]
                rms_ss = rms(BR.T[1])
                '''

                # Pull out data for model
                Vp = Dat_ps['vrel'][ihole_ps,jmid,0,:]
                Vs = Dat_ss['vrel'][ihole_ss,jmid,0,:]
                Pp = Dat_ps['pstat'][ihole_ps,jmid,0,:]
                Ps = Dat_ss['pstat'][ihole_ss,jmid,0,:]
                vpmean = np.tile(np.mean(Vp, axis = 1), (480,1))
                Vp_hat = np.divide(Vp , vpmean.T)
                vsmean = np.tile(np.mean(Vs, axis = 1), (480,1))
                Vs_hat = np.divide(Vs, vsmean.T)
                ppmean = np.tile(np.mean(Pp, axis = 1), (480,1))
                Pp_hat = np.divide(Pp, ppmean.T)
                psmean = np.tile(np.mean(Ps, axis = 1), (480,1))
                Ps_hat = np.divide(Ps, psmean.T)
                x = Dat_ps['x'][:,jmid,0,0]

                #key in form 'Ma_0.70_psi_1.60_phi_0.45'
                #Data['Ma_'+"{:.2f}".format(Mai)+'_psi_'+"{:.2f}".format(Psii)+'_phi_'+"{:.2f}".format(Phii)] = [BR.T[0],BR.T[1],ptp_ps,ptp_ss,rms_ps,rms_ss]
                #find out which way round
                '''
                # Plot the Pressure
                f,a = plt.subplots()  # Create a figure and axis to plot into
                a.plot(ft, Vinf.T)
                a.set_ylabel('Pressure, $Pa$')
                a.set_xlabel('Time, Vane Periods')
                plt.tight_layout()  # Remove extraneous white space
                plt.savefig('Pressure_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai + '.pdf')  # Write out a pdf file

                # Plot the velocity
                f,a = plt.subplots()  # Create a figure and axis to plot into
                a.plot(ft, Vinf.T)
                a.set_ylabel('Velocity, $V$')
                a.set_xlabel('Time, Vane Periods')
                plt.tight_layout()  # Remove extraneous white space
                plt.savefig('Velocity_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai + '.pdf')  # Write out a pdf file               

                # Plot the Blowing ratios
                f,a = plt.subplots()  # Create a figure and axis to plot into
                a.plot(ft, BR.T)
                a.set_ylabel('Hole Blowing Ratio, $BR$')
                a.set_xlabel('Time, Vane Periods')
                plt.tight_layout()  # Remove extraneous white space
                plt.savefig('BR_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai + '.pdf')  # Write out a pdf file
                '''

                f,a = plt.subplots()  # Create a figure and axis to plot into
                for i in range(len(Pp_hat)):
                        position = ihole_ps[i]
                        a.plot(ft,Pp_hat[i],'-', label = 'position = %.2f' %position)  # Plot our data as a new line
                        #a.plot(ft,Vp_hat[1],'--', label = 'position = %.2f' %position)  # Plot our data as a new line
                plt.xlabel('Time, Rotor Periods, $ft$')  # Horizontal axis label
                #plt.ylabel('Static Pressure, $p/\overline{p}$')  # Vertical axis label
                plt.legend()
                plt.tight_layout()  # Remove extraneous white space
                plt.savefig('PS_Velocity_Pressure_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f.pdf' % Mai)  # Write out a pdf file

                f,a = plt.subplots()  # Create a figure and axis to plot into
                for i in range(len(Ps_hat)):
                        position = ihole_ss[i]
                        a.plot(ft,Ps_hat[i],'-', label = 'position = %.2f' %position)  # Plot our data as a new line
                        #a.plot(ft,Vs_hat[i],'--', label = 'position = %.2f' %position)  # Plot our data as a new line
                plt.xlabel('Time, Rotor Periods, $ft$')  # Horizontal axis label
                #plt.ylabel('Static Pressure, $p/\overline{p}$')  # Vertical axis label
                plt.legend()
                plt.tight_layout()  # Remove extraneous white space
                plt.savefig('SS_Velocity_Pressure_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f.pdf' % Mai)  # Write out a pdf file

'''
#looking through phi
#Pressure side peak-to-peak graph
f,a = plt.subplots()
for Phii in [0.45, 0.60, 0.80, 1.00, 1.15]:
        a.plot(Phii, Data['Ma_0.70_psi_1.60_phi_%.2f' %Phii][2], 'x')
a.set_ylabel('Pressure side hole Peak-to-peak Blowing Ratio, $BR$')
a.set_xlabel('Phi')
plt.tight_layout()        
plt.savefig('Pressure_side_peak-to-peak_blowing_ratio_vs_phi.pdf')

#Pressure side rms graph
f,a = plt.subplots()
for Phii in [0.45, 0.60, 0.80, 1.00, 1.15]:
        a.plot(Phii, Data['Ma_0.70_psi_1.60_phi_%.2f' %Phii][4], 'x')
a.set_ylabel('Pressure side hole rms Blowing Ratio, $BR$')
a.set_xlabel('Phi')
plt.tight_layout()        
plt.savefig('Pressure_side_rms_blowing_ratio_vs_phi.pdf')

#Suction side peak-to-peak graph
f,a = plt.subplots()
for Phii in [0.45, 0.60, 0.80, 1.00, 1.15]:
        a.plot(Phii, Data['Ma_0.70_psi_1.60_phi_%.2f' %Phii][3], 'x')
a.set_ylabel('Suction side hole Peak-to-peak Blowing Ratio, $BR$')
a.set_xlabel('Phi')
plt.tight_layout()        
plt.savefig('Suction_side_peak-to-peak_blowing_ratio_vs_phi.pdf')

#Suction side rms graph
f,a = plt.subplots()
for Phii in [0.45, 0.60, 0.80, 1.00, 1.15]:
        a.plot(Phii, Data['Ma_0.70_psi_1.60_phi_%.2f' %Phii][5], 'x')
a.set_ylabel('Suction side hole rms Blowing Ratio, $BR$')
a.set_xlabel('Phi')
plt.tight_layout()        
plt.savefig('Pressure_side_rms_blowing_ratio_vs_phi.pdf')

plt.show()
'''
'''
f,a = plt.subplots()
for x in range(len(Mach)):
        a.plot(ft, BR[x].T, '-', label = 'Mach = %.2f' % Mach[x])
a.set_ylabel('Hole Blowing Ratio, $BR$')
a.set_xlabel('Time, Vane Periods')
plt.tight_layout()  # Remove extraneous white space
plt.savefig('BR_all_Mach.pdf')

plt.show()  # Render the plots
'''       
      