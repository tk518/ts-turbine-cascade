
"""This script reads and plots the blade-to-blade probes."""
import numpy as np  # Multidimensional array library
import probe  # Code for reading TS probe output
import matplotlib.pyplot as plt  # Plotting library
from ts import ts_tstream_reader, ts_tstream_cut  # TS grid reader
from ts import ts_tstream_reader, ts_tstream_patch_kind
import os
from matplotlib.pyplot import cm
#
# Set variables here
#
def rms(x):
        ms = np.sqrt(np.mean(x**2, axis = 0))
        return(ms)

def ptp(x):
    tp = np.max(x, axis = 0) - np.min(x, axis = 0)
    return(tp)

def test_cyclicity(Dat_p, nsteps_cycle, nts):

    absolute_pressure_difference = []
    percentage_pressure_difference = []
    penultimates = nts - nsteps_cycle 
    prepenulatimates = nts - 2*nsteps_cycle

    for point in range(0,nsteps_cycle):
        #compare point in penultimate cycle to final cycle
        #penultimate point - point 1
        point1 = prepenulatimates + point
        #final point - point 2
        point2 = penultimates + point
        difference = Dat_p['pstat'][imid,jmid,0, point1] - Dat_p['pstat'][imid,jmid,0, point2]
        percentage_difference = difference / Dat_p['pstat'][imid,jmid,0, point2]

        #Add results into list
        absolute_pressure_difference.append(difference)
        percentage_pressure_difference.append(abs(percentage_difference))

    #'Average absolute pressure difference = ', sum(absolute_pressure_difference)/len(absolute_pressure_difference)
    #'Average percentage pressure difference = ', sum(percentage_pressure_difference)/len(percentage_pressure_difference)
    #'Maximum percentage cycle difference = ', max(percentage_pressure_difference)*100, '%'
    #'Maximum absolute cycle difference = ', max(absolute_pressure_difference)

    return(absolute_pressure_difference, percentage_pressure_difference)


Phi = [0.80]
Psi = [1.60]
Ma =  [0.70]
Data={}
slip = False

for Psii in Psi:

    for Phii in Phi:

        n = 0
        P_hat1 = []
        P2 = []
        Mach = []

        for Mai in Ma:
            if slip == True:
                output_file_name1 = 'output_2_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f_slip' % Mai  # Location of TS output file
            
            output_file_name = 'output_2_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai

        # We identify a region of the grid using block and patch IDs of stator

        # Load the grid 
        tsr = ts_tstream_reader.TstreamReader()
        g = tsr.read(output_file_name + '.hdf5', read_yplus=True)
    
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
                    if row_str == 'STATOR':
                        probes.append((bid,pid))
                    print('%s, bid=%d, pid=%d, di=%d, dj=%d, dk=%d' % (row_str, bid, pid, di, dj, dk))
        print('***')

        '''
        pid_probe_ps = 9  # Patch ID of probe on pressure side
        pid_probe_ss = 10  # Patch ID of probe on suction side
        '''

        #
        # This next section contains code to read in the data and process it into a
        # convenient form. Only a vague undestanding of this section is needed.
        #

        # Determine number of blades in each row
        '''
        bids = [0,g.get_nb()-1]
        fracann = np.array([g.get_bv('fracann',bi) for bi in bids])
        nblade = np.array([g.get_bv('nblade',bi) for bi in bids])
        nb_row = np.round(fracann * nblade)
        bid_probe = int(nb_row[0])  # Block ID where probes are located
        '''

        # Determine the number of grid points on probe patches
        # (We index the TS grid using i = streamwise, j = spanwise, k = pitchwise)
        p = g.get_patch(probes[0][0],probes[0][1])
        di = p.ien - p.ist
        dj = p.jen - p.jst
        dk = p.ken - p.kst
        probe_shape = [di, dj, dk]  # Numbers of points in i, j, k directions
        print('probe_shape [di, dj, dk]:', probe_shape)

        # Index for the mid-span
        jmid = int(dj/2)

        # Assemble file names for the probes using % substitution
        probe_name = output_file_name + '_probe_%d_%d.dat' % (probes[0][0],probes[0][1])
        #probe_name_ss = output_file_name + '_probe_%d_%d.dat' % (bid_probe,pid_probe_ss)

        # Read the probes
        # The probe data are separate dictionary for each surface of the blade The
        # dictionary is keyed by variable name; the values are numpy arrays with
        # indexes [i = streamwise, j = spanwise, k = pitchwise, n = timewise]
        # For example, to get density at time instant n as a function of
        # axial distance at mid-radius:
        #   Dat_ps['ro'][:,jmid,0,n]
        Dat = probe.read_dat(probe_name, probe_shape)
        #Dat_ss = probe.read_dat(probe_name_ss, probe_shape)

        # Here we extract some parameters from the TS grid to use later
        rpm = g.get_bv('rpm',1)  # RPM in rotor row
        cp = g.get_av('cp')  # Specific heat capacity at const p
        ga = g.get_av('ga')  # Specific heat ratio

        # Get information about time discretisation from TS grid
        freq = g.get_av('frequency')  # Blade passing frequency
        ncycle = g.get_av('ncycle')  # Number of cycles
        nstep_cycle = g.get_av('nstep_cycle')  # Time steps per cycle
        # Individual time step in seconds = blade passing period / steps per cycle
        dt = 1./freq/float(nstep_cycle)
        # Number of time steps = num cycles * steps per cycle
        # nt = ncycle * nstep_cycle
        nt = np.shape(Dat['ro'])[-1]
        print(nt)
        # Make non-dimensional time vector = time in seconds * blade passing frequency
        ft = np.linspace(0.,float(nt-1)*dt,nt) * freq

        # Get secondary vars, things like static pressure, rotor-relative Mach, etc.
        P1, T1 = 1e5, 300.
        Dat = probe.secondary(Dat, rpm, cp, ga, P1, T1)
        #Dat_ss = probe.secondary(Dat_ss, rpm, cp, ga, P1, T1)

        #
        # Finished reading data, now make some plots

        # Streamwise index at mid-chord = number of points in streamwise dirn / 2
        imid = int(di/2)

        # Get static pressure at coordinates of interest
        # i = imid for mid-chord axial location
        # j = jmid for mid-span radial location
        # k = 0 because the patch is at const pitchwise position, on pressure surface
        # n = : for all instants in time 
        P = Dat['pstat'][-1,jmid,:,:]

        # Divide pressure by mean value
        # P is a one-dimensional vector of values of static pressure at each instant in
        # time; np.mean is a function that returns the mean of an array
        P_hat = np.mean(P, axis = 1)


        if slip == True:
                    # Load the grid 
            tsr = ts_tstream_reader.TstreamReader()
            g1 = tsr.read(output_file_name1 + '.hdf5', read_yplus=True)
        
            print('***')
            probes1 = []
            for bid in g.get_block_ids():
                for pid in g.get_patch_ids(bid):
                    patch = g.get_patch(bid,pid)
                    if patch.kind == ts_tstream_patch_kind.probe:
                        rpm = g.get_bv('rpm',bid)
                        row_str = 'STATOR' if rpm==0 else 'ROTOR'
                        di = patch.ien- patch.ist 
                        dj = patch.jen- patch.jst
                        dk = patch.ken- patch.kst
                        if row_str == 'STATOR':
                            probes1.append((bid,pid))
                        print('%s, bid=%d, pid=%d, di=%d, dj=%d, dk=%d' % (row_str, bid, pid, di, dj, dk))
            print('***')

            '''
            pid_probe_ps = 9  # Patch ID of probe on pressure side
            pid_probe_ss = 10  # Patch ID of probe on suction side
            '''

            #
            # This next section contains code to read in the data and process it into a
            # convenient form. Only a vague undestanding of this section is needed.
            #

            # Determine number of blades in each row
            '''
            bids = [0,g.get_nb()-1]
            fracann = np.array([g.get_bv('fracann',bi) for bi in bids])
            nblade = np.array([g.get_bv('nblade',bi) for bi in bids])
            nb_row = np.round(fracann * nblade)
            bid_probe = int(nb_row[0])  # Block ID where probes are located
            '''

            # Determine the number of grid points on probe patches
            # (We index the TS grid using i = streamwise, j = spanwise, k = pitchwise)
            p = g.get_patch(probes[0][0],probes[0][1])
            di1 = p.ien - p.ist
            dj1 = p.jen - p.jst
            dk1 = p.ken - p.kst
            probe_shape1 = [di1, dj1, dk1]  # Numbers of points in i, j, k directions
            print('probe_shape [di, dj, dk]:', probe_shape)

            # Index for the mid-span
            jmid = int(dj1/2)

            # Assemble file names for the probes using % substitution
            probe_name1 = output_file_name1 + '_probe_%d_%d.dat' % (probes1[0][0],probes1[0][1])
            #probe_name_ss = output_file_name + '_probe_%d_%d.dat' % (bid_probe,pid_probe_ss)

            # Read the probes
            # The probe data are separate dictionary for each surface of the blade The
            # dictionary is keyed by variable name; the values are numpy arrays with
            # indexes [i = streamwise, j = spanwise, k = pitchwise, n = timewise]
            # For example, to get density at time instant n as a function of
            # axial distance at mid-radius:
            #   Dat_ps['ro'][:,jmid,0,n]
            Dat1 = probe.read_dat(probe_name1, probe_shape1)
            #Dat_ss = probe.read_dat(probe_name_ss, probe_shape)

            # Here we extract some parameters from the TS grid to use later
            rpm1 = g1.get_bv('rpm',1)  # RPM in rotor row
            cp1 = g1.get_av('cp')  # Specific heat capacity at const p
            ga1 = g1.get_av('ga')  # Specific heat ratio

            # Get information about time discretisation from TS grid
            freq1 = g1.get_av('frequency')  # Blade passing frequency
            ncycle1 = g1.get_av('ncycle')  # Number of cycles
            nstep_cycle1 = g1.get_av('nstep_cycle')  # Time steps per cycle
            # Individual time step in seconds = blade passing period / steps per cycle
            dt1 = 1./freq1/float(nstep_cycle1)
            # Number of time steps = num cycles * steps per cycle
            # nt = ncycle * nstep_cycle
            nt1 = np.shape(Dat1['ro'])[-1]
            print(nt1)
            # Make non-dimensional time vector = time in seconds * blade passing frequency
            ft1 = np.linspace(0.,float(nt1-1)*dt1,nt1) * freq1

            # Get secondary vars, things like static pressure, rotor-relative Mach, etc.
            P1, T1 = 1e5, 300.
            Dat1 = probe.secondary(Dat1, rpm1, cp1, ga1, P1, T1)
            #Dat_ss = probe.secondary(Dat_ss, rpm, cp, ga, P1, T1)

            #
            # Finished reading data, now make some plots

            # Streamwise index at mid-chord = number of points in streamwise dirn / 2
            imid = int(di/2)

            # Get static pressure at coordinates of interest
            # i = imid for mid-chord axial location
            # j = jmid for mid-span radial location
            # k = 0 because the patch is at const pitchwise position, on pressure surface
            # n = : for all instants in time 
            P1 = Dat1['pstat'][-1,jmid,:,:]

            # Divide pressure by mean value
            # P is a one-dimensional vector of values of static pressure at each instant in
            # time; np.mean is a function that returns the mean of an array
            P_hat1 = np.mean(P1, axis = 1)
        # Generate the graph
        '''
        f,a = plt.subplots()  # Create a figure and axis to plot into
        a.plot(ft,P_hat,'-')  # Plot our data as a new line
        plt.xlabel('Time, Rotor Periods, $ft$')  # Horizontal axis label
        plt.ylabel('Static Pressure, $p/\overline{p}$')  # Vertical axis label
        plt.tight_layout()  # Remove extraneous white space
        plt.savefig('unsteady_P_Ma_%.2f.pdf' % Mai)  # Write out a pdf file
        '''
        #
        # Plot time-mean density on pressure side as function of axial location
        #
        path = os.getcwd()
        print(path)
        
        newpath = os.path.join(path, 'psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai)
        try:
            os.mkdir(newpath)
        except OSError:
            print ("Creation of the directory %s failed" % newpath)
        else:
            print ("Successfully created the directory %s " % newpath)

        # Generate the graph
        k = np.linspace(0, dk-1, dk)
        k1 = np.linspace(0, dk1-1, dk1)

        f,a = plt.subplots()  # Create a figure and axis to plot into
        a.plot(P_hat, k,'-', label = 'no slip', '-')  # Plot our data as a new line
        a.plot(P_hat1, k1,'-', label = 'slip', '--')  # Plot our data as a new line
        plt.xlabel('Time-averaged static pressure')  # Horizontal axis label
        #plt.ylabel('Static Pressure, $p/\overline{p}$')  # Vertical axis label
        #plt.legend()
        plt.title('Pressure side plot for Velocity and pressure @ Mach %.2f' % Mai)
        plt.tight_layout()  # Remove extraneous white space
        plt.savefig(newpath + '/Pressure_pitchwise_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f.pdf' % Mai, dpi=200)
        

        # Pull out data for model
        '''
        Vp = Dat_ps['vrel'][:,jmid,0,:]
        Vs = Dat_ss['vrel'][:,jmid,0,:]
        Pp = Dat_ps['pstat'][:,jmid,0,:]
        Ps = Dat_ss['pstat'][:,jmid,0,:]
        Vp_hat = Vp/np.mean(Vp, axis = 0)
        Vs_hat = Vs/np.mean(Vs, axis = 0)
        Pp_hat = Pp/np.mean(Pp, axis = 0)
        Ps_hat = Ps/np.mean(Ps, axis = 0)
        x = Dat_ps['x'][:,jmid,0,0]
        '''

        #Data['Ma_'+"{:.2f}".format(Mai)+'_psi_'+"{:.2f}".format(Psii)+'_phi_'+"{:.2f}".format(Phii)] =[Vp_hat, Vs_hat, Pp_hat, Ps_hat]

        '''
        #For pressure side only to begin
        f,a = plt.subplots()  # Create a figure and axis to plot into

        a.plot(ft,Pp_hat,'-', label = 'Static pressure/Mean static pressure')  # Plot our data as a new line
        a.plot(ft,Vp_hat,'-', label = 'Relative velocity/ Mean relative velocity')  # Plot our data as a new line
        plt.xlabel('Time, Rotor Periods, $ft$')  # Horizontal axis label
        #plt.ylabel('Static Pressure, $p/\overline{p}$')  # Vertical axis label
        plt.legend()
        plt.title('Pressure side plot for Velocity and pressure @ Mach %.2f' % Mai)
        plt.tight_layout()  # Remove extraneous white space
        plt.savefig('Velocity_Pressure_Ma_%.2f.pdf' %Mai)  # Write out a pdf file

        Mach.append(Mai)
        n = n + 1
        '''