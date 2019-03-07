

import time
import glob
from datetime import datetime
from math import floor, log10

import ephem
import numpy as np
from scipy.io import loadmat
import scipy.signal as scisig
import matplotlib.pyplot as plt

from sat_db import active_orbcomm_satellites
from helpers import butter_lowpass_filter, complex_mix


# speed of light
c = 299792458.0 # m/s

# PyEphem observer 
# Imput your receivers latitude, longitude and altitude
# this is the lat/lon that the data was originally recorded at.
lat = 43.802953
lon = -99.210731
alt = 0 
obs = ephem.Observer()
obs.lat, obs.lon = '{}'.format(lat), '{}'.format(lon)
obs.elevation = alt # Technically is the altitude of observer

# Where the data files are located
data_dir = r'./data/'
sample_file = sorted(glob.glob(data_dir + "*.mat"))[0]

# Load the .mat file and print some of the metadata
data = loadmat(sample_file)
print("Timestamp: {}".format(data['timestamp'][0][0]))
print("Satellites in recording: {}".format(', '.join(data['sats'])))
print("Sample rate: {}".format(data['fs'][0][0]))
print("Center frequency: {}".format(data['fc'][0][0]))
frequencies = []
for sat_name in data['sats']:
	freq1, freq2 = active_orbcomm_satellites[sat_name]['frequencies']
	frequencies.append(freq1)
	frequencies.append(freq2)
print("Satellite frequencies: {}".format(', '.join([str(xx) for xx in frequencies])))

# Extract the values for some further processing
samples = data['samples'][0]
center_freq = data['fc'][0][0]
sample_rate = data['fs'][0][0]
timestamp = data['timestamp'][0][0]
obs.date = datetime.utcfromtimestamp(timestamp)


# Plot spectrum of recording
plt.subplot(221)
nperseg = int(sample_rate/100.0)
f, full_pxx = scisig.welch(samples, fs=sample_rate, nperseg=nperseg, scaling='density')
f = (np.roll(f, len(f)/2) + center_freq)/1e6
full_pxx = np.roll(full_pxx, len(full_pxx)/2)
full_pxx = 10*np.log10(full_pxx)
plt.plot(f, full_pxx)
plt.title("Periodogram of recording\nRed dots are orbcomm channels")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
low_point = np.min(full_pxx)
for freq in frequencies:
	plt.plot(freq/1e6, low_point, 'ro')

# Plot one of the channels as baseband
plt.subplot(222)
sat_center_frequency = frequencies[0]

# Use the TLE info that was in the .mat file to calculate doppler shift
# of the satellite's transmission
sat_line0, sat_line1, sat_line2 = [str(xx) for xx in data['tles'][0]]
sat = ephem.readtle(sat_line0, sat_line1, sat_line2)
sat.compute(obs)
relative_vel = sat.range_velocity
doppler = c/(c+relative_vel) * sat_center_frequency - sat_center_frequency

freq_shift = center_freq - sat_center_frequency - doppler
filter_freq = 10e3
mixed_down_samples = complex_mix(samples, freq_shift, sample_rate)
filtered_samples = butter_lowpass_filter(mixed_down_samples, filter_freq, sample_rate, order=5)
f, pxx = scisig.welch(filtered_samples, fs=sample_rate, nperseg=nperseg, scaling='density')
f = (np.roll(f, len(f)/2))
pxx = np.roll(pxx, len(pxx)/2)
pxx = 10*np.log10(pxx)
plt.plot(f, pxx)
plt.xlim([-10e3, 10e3])
plt.ylim([np.min(full_pxx)-1, np.max(pxx)+1])
plt.title("Spectrum of one channel at baseband")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")

# Plot power envelope of the signal
plt.subplot(223)
start = int(1e6)
stop = start + int(100e3)
plt.plot(np.abs(filtered_samples)[start:stop])
plt.title("Power envelope (100k samples)")
plt.xlabel("Sample number")
plt.ylabel("Magnitude")

# Plot complex samples from one channel
plt.subplot(224)
start = int(1e6)
stop = start + int(10e3)
decim = 10
plt.scatter(filtered_samples.real[start:stop:decim], filtered_samples.imag[start:stop:decim])
plt.title("Complex samples (10k samples)")
plt.xlabel("Real")
plt.ylabel("Imag")

plt.tight_layout()
plt.show()