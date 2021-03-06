"""This script reads in an unsteady solution and graphs the results."""
import numpy as np  # Multidimensional array library
import probe  # Code for reading TS probe output
import matplotlib.pyplot as plt  # Plotting library
from ts import ts_tstream_reader  # TS grid reader

#
# Set variables here
#

output_file_name = 'output_2'  # Location of TS output file

# We identify a region of the grid using block and patch IDs
pid_probe_ps = 9  # Patch ID of probe on pressure side
pid_probe_ss = 10  # Patch ID of probe on suction side
bid_probe = 5  # Block ID where probes are located

#
# This next section contains code to read in the data and process it into a
# convenient form. Only a vague undestanding of this section is needed.
#

# Load the grid 
tsr = ts_tstream_reader.TstreamReader()
g = tsr.read(output_file_name + '.hdf5')

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
'''???'''
nt = np.shape(Dat_ps['ro'])[-1] 
print('Number of cycles =', ncycle)
#print('nt =', nt)
# Make non-dimensional time vector = time in seconds * blade passing frequency
ft = np.linspace(0.,float(nt-1)*dt,nt) * freq
#print('ft =', ft)

# Get secondary vars, things like static pressure, rotor-relative Mach, etc.
Dat_ps = probe.secondary(Dat_ps, rpm, cp, ga)
Dat_ss = probe.secondary(Dat_ss, rpm, cp, ga)

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
#print('Length of pressure data =',len(P))
#print('Pressure data =', P)
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
#plt.show()  # Render the plot
plt.savefig('unsteady_P.pdf')  # Write out a pdf file

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
#plt.show()  # Render the plot
plt.savefig('P_x.pdf')  # Write out a pdf file

#
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
#plt.show()  # Render the plot
plt.savefig('unsteady_P.pdf')  # Write out a pdf file


#Point where the change becomes cyclical (number of cycles before the unsteadiness is predictable)
n_cycles_repeating = 50
nt = np.shape(Dat_ps['ro'])[-1] 

penultimate = nt - nstep_cycle 
prepenulatimate = nt - 2*nstep_cycle

P_fourier_penultimate = Dat_ps['pstat'][imid,jmid,0, penultimate:]
P_fourier_prepenultimate = Dat_ps['pstat'][imid,jmid,0, prepenulatimate:penultimate]
#print('length penultimate =', len(P_fourier_penultimate))
#print('length prepenultimate =', len(P_fourier_prepenultimate))
#FFT plot
P_fourier_penultimate = P_fourier_penultimate - np.mean(P_fourier_penultimate)
fourierTransform = np.fft.fft(P_fourier_penultimate)
frequencies = np.fft.fftfreq(len(P_fourier_penultimate) , dt)

#print('Fourier frequencies penultimate = ', frequencies)
#print('fourierTransform penultimate = ', fourierTransform)

#Frequency domain presentation
plt.plot(frequencies, abs(fourierTransform))
plt.title('FFT 79 to 80 passes at Midpoint; Mach = 0.70')
plt.xlabel('Frequency')  
plt.ylabel('Amplitude')
#plt.show()
plt.savefig('FFT_fullCFD_mid.pdf')  # Write out a pdf file

P_fourier_prepenultimate = P_fourier_prepenultimate - np.mean(P_fourier_prepenultimate)
fourierTransform = np.fft.fft(P_fourier_prepenultimate)
frequencies = np.fft.fftfreq(len(P_fourier_prepenultimate) , dt)

#print('Fourier frequencies prepenultimate = ', frequencies)
#print('fourierTransform prepenultimate = ', fourierTransform)

#Frequency domain presentation
plt.plot(frequencies, abs(fourierTransform))
plt.title('FFT 78 to 79 passes at Midpoint; Mach = 0.70')
plt.xlabel('Frequency')  
plt.ylabel('Amplitude')
#plt.show()
plt.savefig('FFT_fullCFD_mid.pdf')  # Write out a pdf file

'''
P2 = P2 - np.mean(P2)
fourierTransform = np.fft.fft(P2)
frequencies = np.fft.fftfreq(len(P2) , dt)

#Frequency domain presentation
plt.plot(frequencies, abs(fourierTransform))
plt.title('FFT 75 to 80 passes at Leading edge; Mach = 0.70 ')
plt.xlabel('Frequency')  
plt.ylabel('Amplitude')
#plt.show()
plt.savefig('FFT_fullCFD_leading.pdf')  # Write out a pdf file

P3 = P3 - np.mean(P3)
fourierTransform = np.fft.fft(P3)
frequencies = np.fft.fftfreq(len(P3) , dt)

#Frequency domain presentation
plt.plot(frequencies, abs(fourierTransform))
plt.title('FFT 75 to 80 passes at Trailing edge; Mach = 0.70')
plt.xlabel('Frequency')  
plt.ylabel('Amplitude')
#plt.show()
plt.savefig('FFT_fullCFD_trailing.pdf')  # Write out a pdf file
'''


#Testing how the level as to which the pressure is repeating.



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

    return(absolute_pressure_difference, percentage_pressure_difference, )




test = test_cyclicity(Dat_ps, nstep_cycle, nt)

print('Average absolute pressure difference = ', sum(test[0])/len(test[0]))
print('Average percentage pressure difference = ', sum(test[1])*100/len(test[1]), '%')
print('Maximum percentage cycle difference = ', max(test[1])*100, '%')
print('Maximum absolute cycle difference = ', max(test[0]))

plt.show()