import numpy as np  # Multidimensional array library
import probe  # Code for reading TS probe output
import matplotlib.pyplot as plt  # Plotting library
#from ts import ts_tstream_reader, ts_tstream_cut  # TS grid reader
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

x = np.linspace(0, 10, 1000)
I = np.sin(x) * np.cos(x[:, np.newaxis])

plt.imshow(I)
plt.colorbar()
plt.cm