##############################################################################
#
# Author: Frank Bieberly
# Date: 30 April 2019
# Name: plot_recording_waterfall.py
# Description:
# This script takes in .mat files (produced by record_spectrum_orbcomm.py) and
# produces a plot of the spectrum over time. This can be useful to find the
# doppler shift of the satellite transmissions. Also, it can help notice
# interference that may be causing problems with reception.
#
##############################################################################

import time
import glob
from math import floor, log10

import numpy as np
from scipy.io import loadmat
import scipy.signal as scisig
import matplotlib.pyplot as plt


data_dir = r'./data/'
sample_files = sorted(glob.glob(data_dir + "*.mat"))

sample_rate = 1.2288e6
nperseg = int(sample_rate/10.0)

waterfall_data = []
for filename in sample_files:

    data = loadmat(filename)
    center_freq = data['fc'][0]
    samples = data['samples'][0]
    f, pxx = scisig.welch(samples, fs=sample_rate, nperseg=nperseg, scaling='density')
    f = (np.roll(f, len(f)/2) + center_freq)/1e6
    pxx = np.roll(pxx, len(pxx)/2)
    pxx = 10*np.log10(pxx)
    waterfall_data.append(pxx)

fig, ax = plt.subplots()
ax.ticklabel_format(useOffset=False, style='plain')
plt.imshow(waterfall_data, aspect='auto',
       animated=True, interpolation=None,extent=[min(f),max(f),0, 100])
# plt.xlim([sat_dict[sat_name]['frequency']/1e6 - 0.02, sat_dict[sat_name]['frequency']/1e6 + 0.02])
plt.colorbar()
plt.title("Spectrum")
plt.ylabel("Time (start at top)")
plt.xlabel("Frequency (MHz)")
plt.show()

