"""This script reads in an unsteady solution and runs the simple hole model."""
import numpy as np  # Multidimensional array library
import probe  # Code for reading TS probe output
import model  # Simple hole model
import matplotlib.pyplot as plt  # Plotting library
from matplotlib.pyplot import cm
from ts import ts_tstream_reader, ts_tstream_patch_kind, ts_tstream_cut  # TS grid reader,TS cutter, TS kind
import os
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


for Psii in Psi:

        for Phii in Phi:
                n = 0
                #Mach = []

                for Mai in Ma:
                        if slip == True:
                            output_file_name = 'output_2_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f_slip' % Mai  # Location of TS output file
                        else:
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

                        print('***')
                        probes = []
                        for bid in g.get_block_ids():
                            for pid in g.get_patch_ids(bid):
                                patch = g.get_patch(bid,pid)
                                if patch.kind == ts_tstream_patch_kind.probe:
                                    rpm = g.get_bv('rpm',bid)
                                    row_str = 'STATOR' if rpm==0 else 'ROTOR'
                                    di = patch.ien- patch.ist 
                                    dj = patch.jen- patch.jst
                                    dk = patch.ken- patch.kst
                                    if row_str == 'ROTOR' and dk != 1:
                                        probes.append((bid,pid, dk))
                                    print('%s, bid=%d, pid=%d, di=%d, dj=%d, dk=%d' % (row_str, bid, pid, di, dj, dk))
                        print('***')

                        print 'probes: ', probes
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

                        pitch_probe = output_file_name + '_probe_%d_%d.dat' % (probes[0][0], probes[0][1])

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
                        rpm = g.get_bv('rpm',probes[0][0])  # RPM in rotor row
                        print 'rpm: ', rpm
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

                        p = g.get_patch(probes[0][0],probes[0][1])
                        di = p.ien - p.ist
                        dj = p.jen - p.jst
                        dk = p.ken - p.kst
                        probe_shape = [di, dj, dk]  # Numbers of points in i, j, k directions

                        Dat_pitch = probe.read_dat(pitch_probe, probe_shape)
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
                        _, P1 = rotor_inlet.mass_avg_1d('pstat')
                        _, Mrel1 = rotor_inlet.mass_avg_1d('mach_rel')
                        _, To1 = rotor_inlet.mass_avg_1d('tstag_rel')
                        _, x1 = rotor_inlet.mass_avg_1d('x')
                        Po1 = model.cf.from_Ma('Po_P',Mrel1,ga)*P1

                        #
                        # Set up the simple hole model
                        #

                        # Choose a constant "pressure margin", or percentage increase of coolant
                        # relative to the inlet stagnation condition
                        PM = 0.04
                        Poc = (1. + PM) * Po1

                        # Fix a stagnation temperature ratio, i.e. the coolant is this much colder than
                        # the main-stream inlet stagnation condition
                        TR = 0.5
                        Toc = TR * To1
                        #print 'Dat_pitch: ', Dat_pitch['rt']
                        jmid_pitchprobe = int(dj/2)
                        pitch_k1 = Dat_pitch['rt'][0,jmid_pitchprobe,-1,0] - Dat_pitch['rt'][0,jmid_pitchprobe,0,0]
                        pitch_k = Dat_ss['rt'][0,jmid,0,0] - Dat_ps['rt'][0,jmid,0,0]
                        print 'pitch length: ', pitch_k
                        print 'pitch length difference: ', pitch_k1-pitch_k
                        x = Dat_ps['x'][:,jmid,0,0]
                        x_hat = (x - x.min())/(x.max() - x.min())
                        U = Dat_ps['U'][0, 0, 0, 0]
                        print 'blade speed: ', Dat_ps['U'][:,0,0,0]
                        #pitch movememt
                        T_offset = pitch_k/U
                        print 'absolute T_offset: ', T_offset
                        # Individual time step in seconds = blade passing period / steps per cycle
                        T_offset = T_offset*freq
                        print 'non-dimensional T_offset: ', T_offset


                        # Pull out data for model
                        # Assume constant Cd
                        Cd = 0.7

                        #pressure side first
                        roinf = Dat_ps_free['ro'][:,jmid,0,:]
                        Vinf = Dat_ps_free['vrel'][:,jmid,0,:]
                        Pinf = Dat_ps_free['pstat'][:,jmid,0,:]
                        # Nondimensionalise data
                        Pinf_Poc, roVinf_Po_cpToc = model.normalise(Poc, Toc, Pinf, roinf, Vinf, cp)
                        # Calculate BR
                        BR_ps = model.evaluate( Pinf_Poc, roVinf_Po_cpToc, Cd, ga )

                        #Suction side
                        roinf = Dat_ss_free['ro'][:,jmid,0,:]
                        Vinf = Dat_ss_free['vrel'][:,jmid,0,:]
                        Pinf = Dat_ss_free['pstat'][:,jmid,0,:]
                        # Nondimensionalise data
                        Pinf_Poc, roVinf_Po_cpToc = model.normalise(Poc, Toc, Pinf, roinf, Vinf, cp)
                        # Calculate BR
                        BR_ss = model.evaluate( Pinf_Poc, roVinf_Po_cpToc, Cd, ga )

                        Pps = Dat_ps_free['pstat'][:,jmid,0,:]
                        Pss = Dat_ss_free['pstat'][:,jmid,0,:]
                        # Finished reading data, now make some plots

                        #key in form 'Ma_0.70_psi_1.60_phi_0.45'
                        if slip == True:
                            Data['Ma_'+"{:.2f}".format(Mai)+'_psi_'+"{:.2f}".format(Psii)+'_phi_'+"{:.2f}".format(Phii)+'_slip'] = [Pps, Pss, BR_ps, BR_ss]
                        else:
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

                        #Make a space-time plot for pressure
                        f,a = plt.subplots()  # Create a figure and axis to plot into
                        pmean_ps = np.divide(Dat_ps_free['pstat'][:,jmid,0,:] - np.tile(np.mean(Dat_ps_free['pstat'][:,jmid,0,:],axis = 1), (480,1)).T, np.tile(np.mean(Dat_ps_free['pstat'][:,jmid,0,:],axis = 1),(480,1)).T)
                        pmean_ss = np.divide(Dat_ss_free['pstat'][:,jmid,0,:] - np.tile(np.mean(Dat_ss_free['pstat'][:,jmid,0,:],axis = 1), (480,1)).T, np.tile(np.mean(Dat_ss_free['pstat'][:,jmid,0,:],axis = 1),(480,1)).T)
                        lev = np.linspace(-0.04,0.04,21)
                        a.contourf(-x_hat, ft - T_offset, pmean_ps.T, lev)
                        a.contourf(x_hat, ft, pmean_ss.T, lev)
                        a.set_ylabel('Time, Period')
                        a.set_xlabel('Suction-Surface   Chord  Pressure-Surface')
                        plt.title('Pressure')
                        plt.colorbar(a.contourf(-x_hat, ft - T_offset, pmean_ps.T, lev))
                        plt.ylim(0,4)
                        plt.tight_layout()
                        if slip == True:
                            plt.savefig(newpath + '/space_time_pressure_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f_slip.pdf' % Mai, dpi=200)
                        else:
                            plt.savefig(newpath + '/space_time_pressure_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f.pdf' % Mai, dpi=200)

                        print 'maximum pressure: ', np.max((np.max(pmean_ps), np.max(pmean_ss)))
                        print 'minimum pressure: ', np.min((np.min(pmean_ps), np.min(pmean_ss)))

                        #Make a space-time plot for velocity
                        f,a = plt.subplots()  # Create a figure and axis to plot into
                        vmean_ps = np.divide(Dat_ps_free['vrel'][:,jmid,0,:] - np.tile(np.mean(Dat_ps_free['vrel'][:,jmid,0,:],axis = 1), (480,1)).T, np.tile(np.mean(Dat_ps_free['vrel'][:,jmid,0,:],axis = 1),(480,1)).T)
                        vmean_ss = np.divide(Dat_ss_free['vrel'][:,jmid,0,:] - np.tile(np.mean(Dat_ss_free['vrel'][:,jmid,0,:],axis = 1), (480,1)).T, np.tile(np.mean(Dat_ss_free['vrel'][:,jmid,0,:],axis = 1),(480,1)).T)
                        lev = np.linspace(-0.12,0.12,21)
                        a.contourf(-x_hat, ft-T_offset, vmean_ps.T, lev)
                        a.contourf(x_hat, ft, vmean_ss.T, lev)
                        a.set_ylabel('Time, Period')
                        a.set_xlabel('Suction-Surface   Chord  Pressure-Surface')
                        plt.title('Velocity')
                        plt.colorbar(a.contourf(-x_hat, ft-T_offset, vmean_ps.T, lev))
                        plt.ylim(0,4)
                        plt.tight_layout()
                        if slip == True:
                            plt.savefig(newpath + '/space_time_velocity_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f_slip.pdf' % Mai, dpi=200)
                        else:
                            plt.savefig(newpath + '/space_time_velocity_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f.pdf' % Mai, dpi=200)

                        print 'maximum velocity: ', np.max((np.max(vmean_ps), np.max(vmean_ss)))
                        print 'minimum velocity: ', np.min((np.min(vmean_ps), np.min(vmean_ss)))
                        
                        #Make a space-time plot for Blowing Ratio
                        f,a = plt.subplots()  # Create a figure and axis to plot into
                        BRmean_ps = np.divide(BR_ps - np.tile(np.mean(BR_ps, axis = 1), (480,1)).T, np.tile(np.mean(BR_ps, axis = 1), (480,1)).T)
                        BRmean_ss = np.divide(BR_ss - np.tile(np.mean(BR_ss, axis = 1), (480,1)).T, np.tile(np.mean(BR_ss, axis = 1), (480,1)).T)
                        lev = np.linspace(-0.1,0.1,21)
                        a.contourf(-x_hat, ft-T_offset, BRmean_ps.T, lev)
                        a.contourf(x_hat, ft, BRmean_ss.T, lev)
                        a.set_ylabel('Time, Period')
                        a.set_xlabel('Suction-Surface   Chord  Pressure-Surface')
                        plt.title('Blowing Ratio')
                        plt.colorbar(a.contourf(-x_hat, ft-T_offset, BRmean_ps.T, lev))
                        plt.ylim(0,4)
                        plt.tight_layout()
                        if slip == True:
                            plt.savefig(newpath + '/space_time_BR_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f_slip.pdf' % Mai, dpi=200)
                        else:
                            plt.savefig(newpath + '/space_time_BR_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f.pdf' % Mai, dpi=200)      

                        print 'maximum BR: ', np.max((np.max(BRmean_ps), np.max(BRmean_ss)))
                        print 'minimum BR: ', np.min((np.min(BRmean_ps), np.min(BRmean_ss))) 