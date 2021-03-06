"""This script reads in an unsteady solution and graphs the results."""
import numpy as np  # Multidimensional array library
import probe  # Code for reading TS probe output
import matplotlib.pyplot as plt  # Plotting library
from ts import ts_tstream_reader  # TS grid reader

#
# Set variables here
#
'''Iterate through all Mach numbers'''

for Mai in [0.70]:

        output_file_name = "output_2"+ '_Ma_%.2f.hdf5' % Mai  # Location of TS output file
        Mach = Mai
        # We identify a region of the grid using block and patch IDs
        pid_probe_ps = 9  # Patch ID of probe on pressure side
        pid_probe_ss = 10  # Patch ID of probe on suction side
        bid_probe = 5  # Block ID where probes are located

        
        # This next section contains code to read in the data and process it into a
        # convenient form. Only a vague undestanding of this section is needed.
        #

        # Load the grid 
       
        tsr = ts_tstream_reader.TstreamReader()
        g = tsr.read(output_file_name)


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

        # Read the probes
        # The probe data are separate dictionary for each surface of the blade The
        # dictionary is keyed by variable name; the values are numpy arrays with
        # indexes [i = streamwise, j = spanwise, k = pitchwise, n = timewise]
        # For example, to get density at time instant n as a function of
        # axial distance at mid-radius:
        #   Dat_ps['ro'][:,jmid,0,n]
        Dat_ps = probe.read_dat(probe_name_ps, probe_shape)
        Dat_ss = probe.read_dat(probe_name_ss, probe_shape)

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
        
        #nt = ncycle * nstep_cycle
        nt = np.shape(Dat_ps['ro'])[-1]
        print(nt)

        # Make non-dimensional time vector = time in seconds * blade passing frequency
        ft = np.linspace(0.,float(nt-1)*dt,nt) * freq

        # Get secondary vars, things like static pressure, rotor-relative Mach, etc.
        Dat_ps = probe.secondary(Dat_ps, rpm, cp, ga)
        Dat_ss = probe.secondary(Dat_ss, rpm, cp, ga)

        '''Finished reading data, now make some plots'''
        
        # Plot static pressure at mid-chord on pressure side edge as function of time

        # Streamwise index at mid-chord = number of points in streamwise dirn / 2
        imid = int(di/2)

        # Get static pressure at coordinates of interest
        # i = imid for mid-chord axial location
        # j = jmid for mid-span radial location
        # k = 0 because the patch is at const pitchwise position, on pressure surface
        # n = : for all instants in time 
        '''Pressure plots at 3 points on a surface vs time'''

        P1 = Dat_ps['pstat'][imid,jmid,0,:]
        P2 = Dat_ps['pstat'][0,jmid,0,:]
        P3 = Dat_ps['pstat'][int(di-1),jmid,0,:]

        # Divide pressure by mean value
        # P is a one-dimensional vector of values of static pressure at each instant in
        # time; np.mean is a function that returns the mean of an array
        P_hat1 = P1 / np.mean(P1)
        P_hat2 = P2 / np.mean(P2)
        P_hat3 = P3 / np.mean(P3)

        # Generate the graph
        f,a = plt.subplots()  # Create a figure and axis to plot into
        a.plot(ft,P_hat1,'-', label = 'Midpoint')  # Plot our data as a new line
        a.plot(ft,P_hat2,'-', label = 'Leading edge')
        a.plot(ft,P_hat3,'-', label = 'Trailing edge')
        plt.legend()

        plt.xlabel('Time, Rotor Periods, $ft$')  # Horizontal axis label
        plt.ylabel('Static Pressure, $p/\overline{p}$')  # Vertical axis label
        plt.tight_layout()  # Remove extraneous white space
        plt.show()  # Render the plot
        plt.savefig('unsteady_P.pdf')  # Write out a pdf file


        '''Plot the gradient plots'''

        f,a = plt.subplots() 
        for Di in [Dat_ps, Dat_ss]:
                P = Di['pstat'][:,jmid,0,:]
                P_hat = P / np.mean(P,axis=1)[0]
                # Get axial coordinates on pressure side 
                # i = : for all axial locations
                # j = jmid for mid-span radial location 
                # k = 0 because the patch is at const pitchwise position, on pressure surface
                # n = 0 for first time step, arbitrary because x is not a function of time.
                x = Di['x'][:,jmid,0,0]
                # Convert to axial chord fraction; use the array min and max functions to get
                # the coordinates at leading and trailing edges respectively.
                x_hat = (x - x.min())/(x.max() - x.min())

                a.plot(x_hat,np.mean(P_hat,axis=1),'-k')  # Plot our data as a new line
                a.plot(x_hat,np.min(P_hat,axis=1),'--k')  # Plot our data as a new line
                a.plot(x_hat,np.max(P_hat,axis=1),'--k')  # Plot our data as a new line

        plt.xlabel('Axial Chord Fraction, $\hat{x}$')  # Horizontal axis label
        # Vertical axis label, start string with r so that \r is not interpreted as a
        # special escape sequence for carriage return
        plt.ylabel(
                r'Static Pressure')
        plt.tight_layout()  # Remove extraneous white space
        plt.show()
        plt.savefig('P_x_Ma_'+str(Mai)+'.pdf')


        # Plot time-mean density on pressure side as function of axial location
        #

        # Get density at coordinates of interest
        # i = : for all axial locations
        # j = jmid for mid-span radial location 
        # k = 0 because the patch is at const pitchwise position, on pressure surface
        # n = : for all instants in time
        ro = Dat_ps['ro'][:,jmid,0,:]

        # Take the time-mean of the density at each axial location
        # ro is a 2D matrix of density values over all axial positions and time steps
        # The first index is i (axial), second index is n (time)
        # We use the np.mean funtion with a keyword argument `axis=1` to specify that we
        # want to take the mean over the second index, i.e. in time and not in space.
        ro_av = np.mean(ro, axis=1)

        # Make non-dimensional with the density at leading edge, at index i=0
        ro_hat = ro_av / ro_av[0]

        # Get axial coordinates on pressure side 
        # i = : for all axial locations
        # j = jmid for mid-span radial location 
        # k = 0 because the patch is at const pitchwise position, on pressure surface
        # n = 0 for first time step, arbitrary because x is not a function of time.
        x = Dat_ps['x'][:,jmid,0,0]

        # Convert to axial chord fraction; use the array min and max functions to get
        # the coordinates at leading and trailing edges respectively.
        x_hat = (x - x.min())/(x.max() - x.min())


        # Other things to try
        #
        #   See https://numpy.org/doc/stable/ for documentation on Numpy
        #
        #   * Frequency spectrum of unsteady pressure at a point, use the Fast
        #   Fourier Transform function np.fft.fft( pressure, axis=?) to get Fourier
        #   coefficients for a series expansion in time. Get frequencies for the bins
        #   using np.fft.fftfreq( len(pressure) , dt)
        #   * Time-mean pressure distribution on pressure and suction sides using
        #   np.mean with the correct axis argument
        #   * Minimum and maximum pressure at each axial location. Use function
        #   np.amax( pressure, axis=? ) to take the maximum value over one index (time)
        #   There is a counterpart np.amin
        #   * Vary the Mach number in `make_design.py` and compare the above for
        #   different Mach numbers


        #Point where the change becomes cyclical (number of cycles before the unsteadiness is predictable)
        n_cycles_repeating = 50

        #FFT plot
        P1 = P1 - np.mean(P1)
        fourierTransform = np.fft.fft(P1)
        frequencies = np.fft.fftfreq(len(P1) , dt)

        #Frequency domain presentation
        plt.plot(frequencies, abs(fourierTransform))
        plt.title('FFT 0 to 80 passes at Midpoint; Mach ='+str(Mach))
        plt.xlabel('Frequency')  
        plt.ylabel('Amplitude')
        plt.show()
        plt.savefig('FFT_fullCFD.pdf')  # Write out a pdf file

        P2 = P2 - np.mean(P2)
        fourierTransform = np.fft.fft(P2)
        frequencies = np.fft.fftfreq(len(P2) , dt)

        #Frequency domain presentation
        plt.plot(frequencies, abs(fourierTransform))
        plt.title('FFT 0 to 80 passes at Leading edge; Mach ='+str(Mach))
        plt.xlabel('Frequency')  
        plt.ylabel('Amplitude')
        plt.show()
        plt.savefig('FFT_fullCFD.pdf')  # Write out a pdf file

        P3 = P3 - np.mean(P3)
        fourierTransform = np.fft.fft(P3)
        frequencies = np.fft.fftfreq(len(P3) , dt)

        #Frequency domain presentation
        plt.plot(frequencies, abs(fourierTransform))
        plt.title('FFT 0 to 80 passes at Trailing edge; Mach ='+str(Mach))
        plt.xlabel('Frequency')  
        plt.ylabel('Amplitude')
        plt.show()
        plt.savefig('FFT_fullCFD.pdf')  # Write out a pdf file

        #Time averaged pressure vs maximum and minimum pressure
        # Get pressure at coordinates of interest
        # i = : for all axial locations
        # j = jmid for mid-span radial location 
        # k = 0 because the patch is at const pitchwise position, on pressure surface
        # n = : for all instants in time
        Pps = Dat_ps['pstat'][:,jmid,0,:]
        Pss = Dat_ss['pstat'][:,jmid,0,:]

        # Take the time-mean of the pressure at each axial location
        # P is a 2D matrix of density values over all axial positions and time steps
        # The first index is i (axial), second index is n (time)
        # We use the np.mean funtion with a keyword argument `axis=1` to specify that we
        # want to take the mean over the second index, i.e. in time and not in space.
        Pps_av = np.mean(Pps, axis=1)
        Ppsmax=[]
        for array in Pps:
                Ppsmax.append(max(array))
        Ppsmin=[]
        for array in Pps:
                Ppsmin.append(min(array))
        # Make non-dimensional with the pressure at leading edge, at index i=0
        Pps_hat = Pps_av / Pps_av[0]
        Ppsmax_hat = Ppsmax / Ppsmax[0]
        Ppsmin_hat = Ppsmin / Ppsmin[0]

        Pss_av = np.mean(Pss, axis=1)
        Pssmax=[]
        for array in Pss:
                Pssmax.append(max(array))
        Pssmin=[]
        for array in Pss:
                Pssmin.append(min(array))
        # Make non-dimensional with the pressure at leading edge, at index i=0
        Pss_hat = Pss_av / Pss_av[0]
        Pssmax_hat = Pssmax / Pssmax[0]
        Pssmin_hat = Pssmin / Pssmin[0]
        # Get axial coordinates on pressure side 
        # i = : for all axial locations
        # j = jmid for mid-span radial location 
        # k = 0 because the patch is at const pitchwise position, on pressure surface
        # n = 0 for first time step, arbitrary because x is not a function of time.
        xps = Dat_ps['x'][:,jmid,0,0]
        xss = Dat_ss['x'][:,jmid,0,0]

        # Convert to axial chord fraction; use the array min and max functions to get
        # the coordinates at leading and trailing edges respectively.
        xps_hat = (xps - xps.min())/(xps.max() - xps.min())
        xss_hat = (xss - xss.min())/(xss.max() - xss.min())

        f,a = plt.subplots()  # Create a figure and axis to plot into

        a.plot(xps_hat,Pps_hat,'-', label = 'average pressure - pressure side')  # Plot our data as a new line
        a.plot(xps_hat,Ppsmax_hat, '-',label = 'max pressure - pressure side')
        a.plot(xps_hat,Ppsmin_hat, '-', label = 'min pressure - pressure side')

        a.plot(xss_hat,Pss_hat,'-', label = 'average pressure - suction side')  # Plot our data as a new line
        a.plot(xss_hat,Pssmax_hat, '-',label = 'max pressure - suction side')
        a.plot(xss_hat,Pssmin_hat, '-', label = 'min pressure - suction side')

        plt.legend()

        plt.title('Max and Min pressures from 0 to 80 passes; Mach ='+str(Mach))

        plt.xlabel('Axial Chord Fraction, $\hat{x}$')  # Horizontal axis label
        # Vertical axis label, start string with r so that \r is not interpreted as a
        # special escape sequence for carriage return
        plt.ylabel(
                r'Time-averaged Static Pressure, $p/\overline{p}$')
        plt.tight_layout()  # Remove extraneous white space
        plt.show()  # Render the plot
        plt.savefig('P_x.pdf')  # Write out a pdf file


        Pps = Dat_ps['pstat'][:,jmid,0,nstep_cycle*n_cycles_repeating:]
        Pss = Dat_ss['pstat'][:,jmid,0,nstep_cycle*n_cycles_repeating:]

        # Take the time-mean of the pressure at each axial location
        # P is a 2D matrix of density values over all axial positions and time steps
        # The first index is i (axial), second index is n (time)
        # We use the np.mean funtion with a keyword argument `axis=1` to specify that we
        # want to take the mean over the second index, i.e. in time and not in space.
        Pps_av = np.mean(Pps, axis=1)
        Ppsmax=[]
        for array in Pps:
                Ppsmax.append(max(array))
        Ppsmin=[]
        for array in Pps:
                Ppsmin.append(min(array))
        # Make non-dimensional with the pressure at leading edge, at index i=0
        Pps_hat = Pps_av / Pps_av[0]
        Ppsmax_hat = Ppsmax / Ppsmax[0]
        Ppsmin_hat = Ppsmin / Ppsmin[0]

        Pss_av = np.mean(Pss, axis=1)
        Pssmax=[]
        for array in Pss:
                Pssmax.append(max(array))
        Pssmin=[]
        for array in Pss:
                Pssmin.append(min(array))
        # Make non-dimensional with the pressure at leading edge, at index i=0
        Pss_hat = Pss_av / Pss_av[0]
        Pssmax_hat = Pssmax / Pssmax[0]
        Pssmin_hat = Pssmin / Pssmin[0]
        # Get axial coordinates on pressure side 
        # i = : for all axial locations
        # j = jmid for mid-span radial location 
        # k = 0 because the patch is at const pitchwise position, on pressure surface
        # n = 0 for first time step, arbitrary because x is not a function of time.
        xps = Dat_ps['x'][:,jmid,0,0]
        xss = Dat_ss['x'][:,jmid,0,0]

        # Convert to axial chord fraction; use the array min and max functions to get
        # the coordinates at leading and trailing edges respectively.
        xps_hat = (xps - xps.min())/(xps.max() - xps.min())
        xss_hat = (xss - xss.min())/(xss.max() - xss.min())

        f,a = plt.subplots()  # Create a figure and axis to plot into

        a.plot(xps_hat,Pps_hat,'-', label = 'average pressure - pressure side')  # Plot our data as a new line
        a.plot(xps_hat,Ppsmax_hat, '-',label = 'max pressure - pressure side')
        a.plot(xps_hat,Ppsmin_hat, '-', label = 'min pressure - pressure side')

        a.plot(xss_hat,Pss_hat,'-', label = 'average pressure - suction side')  # Plot our data as a new line
        a.plot(xss_hat,Pssmax_hat, '-',label = 'max pressure - suction side')
        a.plot(xss_hat,Pssmin_hat, '-', label = 'min pressure - suction side')

        plt.legend()

        plt.title('Max and Min pressures from ' +str(n_cycles_repeating)+ ' to 80 passes; Mach =' +str(Mach))

        plt.xlabel('Axial Chord Fraction, $\hat{x}$')  # Horizontal axis label
        # Vertical axis label, start string with r so that \r is not interpreted as a
        # special escape sequence for carriage return
        plt.ylabel(
                r'Time-averaged Static Pressure, $p/\overline{p}$')
        plt.tight_layout()  # Remove extraneous white space
        plt.show()  # Render the plot
        plt.savefig('P_x_cyclicalpart.pdf')  # Write out a pdf file

        P1 = Dat_ps['pstat'][imid,jmid,0,nstep_cycle*n_cycles_repeating::]
        P2 = Dat_ps['pstat'][0,jmid,0,nstep_cycle*n_cycles_repeating::]
        P3 = Dat_ps['pstat'][int(di-1),jmid,0,nstep_cycle*n_cycles_repeating::]

        # Divide pressure by mean value
        # P is a one-dimensional vector of values of static pressure at each instant in
        # time; np.mean is a function that returns the mean of an array
        P_hat1 = P1 / np.mean(P1)
        P_hat2 = P2 / np.mean(P2)
        P_hat3 = P3 / np.mean(P3)
        # Make non-dimensional time vector = time in seconds * blade passing frequency
        ft1 = np.linspace(float(nstep_cycle*n_cycles_repeating)*dt,float(nt-1)*dt,nt - nstep_cycle*n_cycles_repeating) * freq

        # Generate the graph
        f,a = plt.subplots()  # Create a figure and axis to plot into
        a.plot(ft1,P_hat1,'-', label = 'Midpoint')  # Plot our data as a new line
        a.plot(ft1,P_hat2,'-', label = 'Leading edge')
        a.plot(ft1,P_hat3,'-', label = 'Trailing edge')
        plt.legend()
        plt.title('pressure variation at 3 points: from' +str(n_cycles_repeating)+ ' to 80 passes; Mach ='+str(Mach))

        plt.xlabel('Time, Rotor Periods, $ft$')  # Horizontal axis label
        plt.ylabel('Static Pressure, $p/\overline{p}$')  # Vertical axis label
        plt.tight_layout()  # Remove extraneous white space
        plt.show()  # Render the plot
        plt.savefig('unsteady_P_cyclicalpart.pdf')  # Write out a pdf file

        P1 = P1 - np.mean(P1)
        fourierTransform = np.fft.fft(P1)
        frequencies = np.fft.fftfreq(len(P1) , dt)

        print('Fourier frequencies = ', frequencies)
        print('fourierTransform = ', fourierTransform)

        #Frequency domain presentation
        plt.plot(frequencies, abs(fourierTransform))
        plt.title('FFT ' +str(n_cycles_repeating)+ ' to 80 passes at midpoint; Mach ='+str(Mach))
        plt.xlabel('Frequency')  
        plt.ylabel('Amplitude')
        plt.show()
        plt.savefig('FFT_cyclicalCFD.pdf')

        P2 = P2 - np.mean(P2)
        fourierTransform = np.fft.fft(P2)
        frequencies = np.fft.fftfreq(len(P2) , dt)

'''
        #Frequency domain presentation
        plt.plot(frequencies, abs(fourierTransform))
        plt.title('FFT ' +str(n_cycles_repeating)+ ' to 80 passes at Leading edge; Mach ='+str(Mach))
        plt.xlabel('Frequency')  
        plt.ylabel('Amplitude')
        plt.show()
        plt.savefig('FFT_fullCFD.pdf')  # Write out a pdf file

        P3 = P3 - np.mean(P3)
        fourierTransform = np.fft.fft(P3)
        frequencies = np.fft.fftfreq(len(P3) , dt)

        #Frequency domain presentation
        plt.plot(frequencies, abs(fourierTransform))
        plt.title('FFT ' +str(n_cycles_repeating)+ ' 80 passes at Trailing edge; Mach ='+str(Mach))
        plt.xlabel('Frequency')  
        plt.ylabel('Amplitude')
        plt.show()
        plt.savefig('FFT_fullCFD.pdf')  # Write out a pdf file

'''
