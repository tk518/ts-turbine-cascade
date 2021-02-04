"""This script reads in an unsteady solution and graphs the results."""
import numpy as np  # Multidimensional array library
import probe  # Code for reading TS probe output
import matplotlib.pyplot as plt  # Plotting library
from ts import ts_tstream_reader  # TS grid reader

#
# Set variables here
#

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


Phi = [0.40, 0.60, 0.80, 1.00, 1.20]
Psi = [0.80, 1.20, 1.60, 2.00, 2.40]
Ma =  [0.70]

for Psii in Psi:

    for Phii in Phi:

        n = 0
        P_hat1 = []
        P1 = []
        Mach = []

        for Mai in Ma:

            output_file_name = 'output_2_psi_%.2f' %Psii + '_phi_%.2f' %Phii + '_Ma_%.2f' % Mai  # Location of TS output file

            # We identify a region of the grid using block and patch IDs
            pid_probe_ps = 9  # Patch ID of probe on pressure side
            pid_probe_ss = 10  # Patch ID of probe on suction side

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
            # nt = ncycle * nstep_cycle
            nt = np.shape(Dat_ps['ro'])[-1]
            print(nt)
            # Make non-dimensional time vector = time in seconds * blade passing frequency
            ft = np.linspace(0.,float(nt-1)*dt,nt) * freq

            # Get secondary vars, things like static pressure, rotor-relative Mach, etc.
            P1, T1 = 1e5, 300.
            Dat_ps = probe.secondary(Dat_ps, rpm, cp, ga, P1, T1)
            Dat_ss = probe.secondary(Dat_ss, rpm, cp, ga, P1, T1)

            #
            # Finished reading data, now make some plots
            #

            #
            # Plot static pressure at mid-chord on pressure side edge as function of time
            #

            # Streamwise index at mid-chord = number of points in streamwise dirn / 2
            imid = int(di/2)

            # Get static pressure at coordinates of interest
            # i = imid for mid-chord axial location
            # j = jmid for mid-span radial location
            # k = 0 because the patch is at const pitchwise position, on pressure surface
            # n = : for all instants in time 
            P = Dat_ps['pstat'][0,jmid,0,:]

            # Divide pressure by mean value
            # P is a one-dimensional vector of values of static pressure at each instant in
            # time; np.mean is a function that returns the mean of an array
            P_hat = P / np.mean(P)

            # Generate the graph
            f,a = plt.subplots()  # Create a figure and axis to plot into
            a.plot(ft,P_hat,'-')  # Plot our data as a new line
            plt.xlabel('Time, Rotor Periods, $ft$')  # Horizontal axis label
            plt.ylabel('Static Pressure, $p/\overline{p}$')  # Vertical axis label
            plt.tight_layout()  # Remove extraneous white space
            plt.savefig('unsteady_P_Ma_%.2f.pdf' % Mai)  # Write out a pdf file

            #
            # Plot time-mean density on pressure side as function of axial location
            #

            # Generate the graph
            f,a = plt.subplots()  # Create a figure and axis to plot into

            # Get density at coordinates of interest
            # i = : for all axial locations
            # j = jmid for mid-span radial location 
            # k = 0 because the patch is at const pitchwise position, on pressure surface
            # n = : for all instants in time
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
            plt.savefig('P_x_Ma_%.2f.pdf' % Mai)  # Write out a pdf file

            P1.append(Dat_ps['pstat'][imid,jmid,0,:])
            P_hat1.append(P1[n] / np.mean(P1[n]))

            # Choose a hole position
            ihole_ps = 20
            ihole_ss = 40

            # Pull out data for model

            Vp = Dat_ps['vrel'][ihole_ps,jmid,0,:]
            Vs = Dat_ss['vrel'][ihole_ss,jmid,0,:]
            Pp = Dat_ps['pstat'][ihole_ps,jmid,0,:]
            Ps = Dat_ss['pstat'][ihole_ss,jmid,0,:]
            Vp_hat = Vp/np.mean(Vp)
            Vs_hat = Vs/np.mean(Vs)
            Pp_hat = Pp/np.mean(Pp)
            Ps_hat = Ps/np.mean(Ps)

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

            test = test_cyclicity(Dat_ps, nstep_cycle, nt)

            print('Average absolute pressure difference = at Mach %.2f' % Mai, sum(test[0])/len(test[0]))
            print('Average percentage pressure difference = at Mach %.2f' % Mai, sum(test[1])*100/len(test[1]), '%')
            print('Maximum percentage cycle difference = at Mach %.2f' % Mai, max(test[1])*100, '%')
            print('Maximum absolute cycle difference = at Mach %.2f' % Mai, max(test[0]))


plt.show()  # Render the plots

f,a = plt.subplots()  # Create a figure and axis to plot into

for x in range(0, n):
    a.plot(ft,P_hat1[x],'-', label = 'Midpoint @ Mach = %.2f' % Mach[x])  # Plot our data as a new line

plt.xlabel('Time, Rotor Periods, $ft$')  # Horizontal axis label
plt.ylabel('Static Pressure, $p/\overline{p}$')  # Vertical axis label
plt.legend()
plt.tight_layout()  # Remove extraneous white space
plt.savefig('unsteady_P_at_different_Mach_No.pdf')  # Write out a pdf file
plt.show()

#I want to see 
#a) blowing ratio for multiple Mach numbers [x]
#b) Static pressure change and velocity at the same point on the same graph [x]
