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

Phi = [0.40, 0.60, 0.80, 1.00, 1.20]
Psi = [0.80, 1.20, 1.60, 2.00, 2.40]
Ma =  [0.70]

slip = True

if slip == True:
        for Psii in Psi:

                for Phii in Phi:
                        n = 0
                        #Mach = []

                        for Mai in Ma:

                                output_file_name = 'output_2_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f_slip' % Mai  # Location of TS output file

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
                                print(x.min())
                                print(x.max())

                                # Choose a hole position
                                ihole_ps = x
                                ihole_ss = 40

                                # Pull out data for model

                                roinf = np.stack((Dat_ps_free['ro'][:,jmid,0,:],
                                                Dat_ss_free['ro'][:,jmid,0,:]))
                                Vinf = np.stack((Dat_ps_free['vrel'][:,jmid,0,:],
                                                Dat_ss_free['vrel'][:,jmid,0,:]))
                                Pinf = np.stack((Dat_ps_free['pstat'][:,jmid,0,:],
                                                Dat_ss_free['pstat'][:,jmid,0,:]))

                                # Nondimensionalise data
                                Pinf_Poc, roVinf_Po_cpToc = model.normalise(Poc, Toc, Pinf, roinf, Vinf, cp)

                                # Assume constant Cd
                                Cd = 0.7

                                # Calculate BR
                                BR = model.evaluate( Pinf_Poc, roVinf_Po_cpToc, Cd, ga )
                                #print 'BR: ', BR
                                #
                                # Finished reading data, now make some plots
                                #
                                '''
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
                                        plt.savefig('hole_posn.pdf')  # Write out a pdf file

                                '''
                                #peak to peak amplitude

                                #pressure side, if pressure side is [0]
                                ptp_ps = np.max(BR[0].T, axis = 0) - np.min(BR[0].T, axis = 0)
                                #suction side, if pressure side is [1]
                                ptp_ss = np.max(BR[1].T, axis = 0) - np.min(BR[1].T, axis = 0)

                                #rms BR
                                #pressure side, if pressure side is [0]
                                rms_ps = rms(BR[0].T)
                                #suction side, if pressure side is [1]
                                rms_ss = rms(BR[1].T)

                                #key in form 'Ma_0.70_psi_1.60_phi_0.45'
                                Data['Ma_'+"{:.2f}".format(Mai)+'_psi_'+"{:.2f}".format(Psii)+'_phi_'+"{:.2f}".format(Phii)+'_slip'] = [BR[0].T,BR[1].T,ptp_ps,ptp_ss,rms_ps,rms_ss]
                                #find out which way round

                                #print 'ptp_ps: ', ptp_ps
                                print 'x: ', x
                                print 'rns_ps: ',rms_ps
                                print 'len x: ', len(x), 'len rms: ', len(rms_ps)		

                                # Plot the peak-to-peak Blowing ratios on pressure surface along chord
                                f,a = plt.subplots()  # Create a figure and axis to plot into
                                a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f_slip' %Phii][2])
                                a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f_slip' %Phii][3])
                                a.set_ylabel('Peak-to-Peak Blowing Ratio, $BR$')
                                a.set_xlabel('Chord, x')
                                plt.tight_layout()  # Remove extraneous white space
                                plt.savefig('BR_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai + '_along_chord_slip.pdf')  # Write out a pdf file

                                n = n + 1
                
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
                        print(x.min())
                        print(x.max())

                        # Choose a hole position
                        ihole_ps = x
                        ihole_ss = 40

                        # Pull out data for model

                        roinf = np.stack((Dat_ps_free['ro'][:,jmid,0,:],
                                        Dat_ss_free['ro'][:,jmid,0,:]))
                        Vinf = np.stack((Dat_ps_free['vrel'][:,jmid,0,:],
                                        Dat_ss_free['vrel'][:,jmid,0,:]))
                        Pinf = np.stack((Dat_ps_free['pstat'][:,jmid,0,:],
                                        Dat_ss_free['pstat'][:,jmid,0,:]))

                        # Nondimensionalise data
                        Pinf_Poc, roVinf_Po_cpToc = model.normalise(Poc, Toc, Pinf, roinf, Vinf, cp)

                        # Assume constant Cd
                        Cd = 0.7

                        # Calculate BR
                        BR = model.evaluate( Pinf_Poc, roVinf_Po_cpToc, Cd, ga )
                        #print 'BR: ', BR
                        #
                        # Finished reading data, now make some plots
                        #
                        '''
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
                                plt.savefig('hole_posn.pdf')  # Write out a pdf file

                        '''
                        #peak to peak amplitude

                        #pressure side, if pressure side is [0]
                        ptp_ps = np.max(BR[0].T, axis = 0) - np.min(BR[0].T, axis = 0)
                        #suction side, if pressure side is [1]
                        ptp_ss = np.max(BR[1].T, axis = 0) - np.min(BR[1].T, axis = 0)

                        #rms BR
                        #pressure side, if pressure side is [0]
                        rms_ps = rms(BR[0].T)
                        #suction side, if pressure side is [1]
                        rms_ss = rms(BR[1].T)

                        #key in form 'Ma_0.70_psi_1.60_phi_0.45'
                        Data['Ma_'+"{:.2f}".format(Mai)+'_psi_'+"{:.2f}".format(Psii)+'_phi_'+"{:.2f}".format(Phii)] = [BR[0].T,BR[1].T,ptp_ps,ptp_ss,rms_ps,rms_ss]
                        #find out which way round

                        #print 'ptp_ps: ', ptp_ps
                        print 'x: ', x
                        print 'rns_ps: ',rms_ps
                        print 'len x: ', len(x), 'len rms: ', len(rms_ps)		

                        # Plot the peak-to-peak Blowing ratios on pressure surface along chord
                        f,a = plt.subplots()  # Create a figure and axis to plot into
                        a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f' %Phii][2])
			a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f' %Phii][3])
                        a.set_ylabel('Peak-to-Peak Blowing Ratio, $BR$')
                        a.set_xlabel('Chord, x')
                        plt.tight_layout()  # Remove extraneous white space
                        plt.savefig('BR_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai + '_along_chord.pdf')  # Write out a pdf file

                        n = n + 1


# looking through phi
# Pressure side peak-to-peak graph

# print 'difference between slip and normal ps: ', Data['Ma_0.70_psi_1.60_phi_0.45'][2] - Data['Ma_0.70_psi_1.60_phi_0.45_slip'][2]
# print 'difference between slip and normal ss: ', Data['Ma_0.70_psi_1.60_phi_0.45'][3] - Data['Ma_0.70_psi_1.60_phi_0.45_slip'][3]

f,a = plt.subplots()

n = len(Phi)
color=iter(cm.rainbow(np.linspace(0,1,n)))
for Phii in Phi:
	c = next(color)
	a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f' %Phii][2], '-', label = 'BR @ Phi = %.2f' %Phii, c=c)
        if slip == True:
                a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f_slip' %Phii][2], '--', label = 'BR @ Phi = %.2f with slip' %Phii, c=c)
plt.ylabel('Pressure side Peak-to-peak Blowing Ratio, $BR$')
plt.xlabel('Axial displacement, $x$')
plt.tight_layout()
#plt.ylim(0,0.9)   
plt.legend(loc="best", ncol=2)
if slip == True:
        plt.savefig('Pressure_side_peak-to-peak_blowing_ratio_vs_x_slip.pdf')
else:
        plt.savefig('Pressure_side_peak-to-peak_blowing_ratio_vs_x.pdf')

#Pressure side rms graph
f,a = plt.subplots()
n = len(Phi)
color=iter(cm.rainbow(np.linspace(0,1,n)))
for Phii in Phi:
	c = next(color)
        a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f' %Phii][4], '-', label = 'BR @ Phi = %.2f' %Phii, c=c)
        if slip == True:
                a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f_slip' %Phii][2], '--', label = 'BR @ Phi = %.2f with slip' %Phii, c=c)
a.set_ylabel('Pressure side hole rms Blowing Ratio, $BR$')
a.set_xlabel('Axial displacement, $x$')
plt.tight_layout()  
plt.legend(loc="best", ncol=2)     
if slip == True: 
        plt.savefig('Pressure_side_rms_blowing_ratio_vs_x_slip.pdf')
else:
        plt.savefig('Pressure_side_rms_blowing_ratio_vs_x.pdf')

#Suction side peak-to-peak graph
f,a = plt.subplots()
n = len(Phi)
color=iter(cm.rainbow(np.linspace(0,1,n)))
for Phii in Phi:
	c = next(color)
        a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f' %Phii][3], '-', label = 'BR @ Phi = %.2f' %Phii, c=c)
        if slip == True:
                a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f_slip' %Phii][2], '--', label = 'BR @ Phi = %.2f with slip' %Phii, c=c)
plt.ylabel('Suction side Peak-to-peak Blowing Ratio, $BR$')
plt.xlabel('Chord, $x$')
plt.tight_layout()
#plt.ylim(0,0.9)   
plt.legend(loc="best", ncol=2)
if slip == True:
        plt.savefig('Suction_side_peak-to-peak_blowing_ratio_vs_x_slip.pdf')
else:
        plt.savefig('Suction_side_peak-to-peak_blowing_ratio_vs_x.pdf')

#Suction side rms graph
f,a = plt.subplots()
n = len(Phi)
color=iter(cm.rainbow(np.linspace(0,1,n)))
for Phii in Phi:
	c = next(color)
        a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f' %Phii][5], '-', label = 'BR @ Phi = %.2f' %Phii, c=c)
        if slip == True:
                a.plot(x, Data['Ma_0.70_psi_1.60_phi_%.2f_slip' %Phii][2], '--', label = 'BR @ Phi = %.2f with slip' %Phii, c=c)
a.set_ylabel('Suction side rms Blowing Ratio, $BR$')
a.set_xlabel('Axial displacement, $x$')
plt.tight_layout()        
plt.legend(loc="best", ncol=2)
if slip == True:
        plt.savefig('Suction_side_rms_blowing_ratio_vs_x_slip.pdf')
else:
        plt.savefig('Suction_side_rms_blowing_ratio_vs_x.pdf')

plt.show()

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