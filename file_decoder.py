

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

# Get the TLE information from the .mat file
sat_line0, sat_line1, sat_line2 = [str(xx) for xx in data['tles'][0]]
sat = ephem.readtle(sat_line0, sat_line1, sat_line2)
sat.compute(obs)

# Decode the lower channel
sat_center_frequency = frequencies[0]

# Use the TLE info that was in the .mat file to calculate doppler shift
# of the satellite's transmission
relative_vel = sat.range_velocity
doppler = c/(c+relative_vel) * sat_center_frequency - sat_center_frequency

# Mix the samples to baseband (compensating for doppler shift)
# There will be a residual carrier error because of the RLTSDR frequency offset
freq_shift = center_freq - sat_center_frequency - doppler
mixed_down_samples = complex_mix(samples, freq_shift, sample_rate)

# Low pass filter the signal
filter_freq = 10e3
filtered_samples = butter_lowpass_filter(mixed_down_samples, filter_freq, sample_rate, order=5)

# Decimated signal
baud_rate = 4800.0
samples_per_symbol = 8
decimation = int(sample_rate/(samples_per_symbol*baud_rate))
decimated_samples = filtered_samples[::decimation]

# estimate remaining carrier error (RTLSDR frequency error)
# signal to the fourth power, take the fft, peak is at frequency offset
rbw = 10
nperseg = int(sample_rate/decimation/rbw)
signal_to_4th_power = np.power(decimated_samples, 4)
f, pxx = scisig.welch(signal_to_4th_power, fs=sample_rate/decimation, nperseg=nperseg, scaling='density')
f = (np.roll(f, len(f)/2))
pxx = np.roll(pxx, len(pxx)/2)
search_window = int(1000 / ((sample_rate/decimation)/nperseg))	# search +/- 1 kHz around fc
frequency_peak = np.argmax(pxx[len(pxx)/2 - search_window:len(pxx)/2 + search_window])
freq_offset = -(frequency_peak - search_window)*(sample_rate/decimation/nperseg) / 4
baseband_samples = complex_mix(decimated_samples, freq_offset, sample_rate/decimation)


# Costas loop carrier/phase recovery
filt_order = 50	# number of taps (plus 1)
filt_bands = np.array([0.0, 0.05, 0.10, 1.0]) * (sample_rate/decimation)/2
filt_amp = np.array([1.0, 0.0])
fir_taps = np.flip(scisig.remez(filt_order + 1, filt_bands, filt_amp, fs=sample_rate/decimation), 0)

# freq, response = scisig.freqz(fir_taps)
# plt.figure()
# plt.semilogy(0.5*sample_rate/decimation*freq/np.pi, np.abs(response), 'b-')
# plt.grid(alpha=0.25)
# plt.title("Frequency response of LPF in Costas loop")
# plt.xlabel('Frequency (Hz)')
# plt.ylabel('Gain')

mu = 0.005	# update gain
f0 = 0		# expected baseband frequency
phase_est = np.zeros(len(baseband_samples)+1)

z1 = np.zeros(filt_order + 1, dtype=np.complex64)
z2 = np.zeros(filt_order + 1, dtype=np.complex64)
z3 = np.zeros(filt_order + 1, dtype=np.complex64)
z4 = np.zeros(filt_order + 1, dtype=np.complex64)

for idx, sample in enumerate(baseband_samples):
	sample = 2*sample

	z1[:-1] = z1[1:]
	z1[-1] = sample * np.exp(1j * 2 * np.pi * f0 * idx * 1/(sample_rate/decimation) + phase_est[idx])
	z2[:-1] = z2[1:]
	z2[-1] = sample * np.exp(1j * 2 * np.pi * f0 * idx * 1/(sample_rate/decimation) + phase_est[idx] + np.pi/4)
	z3[:-1] = z3[1:]
	z3[-1] = sample * np.exp(1j * 2 * np.pi * f0 * idx * 1/(sample_rate/decimation) + phase_est[idx] + np.pi/2)
	z4[:-1] = z4[1:]
	z4[-1] = sample * np.exp(1j * 2 * np.pi * f0 * idx * 1/(sample_rate/decimation) + phase_est[idx] + np.pi*3/4)

	lpf1 = np.dot(fir_taps, z1)
	lpf2 = np.dot(fir_taps, z2)
	lpf3 = np.dot(fir_taps, z3)
	lpf4 = np.dot(fir_taps, z4)

	phase_est[idx+1] = phase_est[idx] + mu * np.angle(lpf1 * lpf2 * lpf3 * lpf4)

plt.figure()
plt.plot(phase_est)
plt.title("Costas Loop Phase estimate")
plt.xlabel("Sample number")
plt.ylabel("Phase")


# Plot spectrum of recording
plt.figure()
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
#plot the decimated signal
f, pxx = scisig.welch(baseband_samples, fs=sample_rate/decimation, nperseg=nperseg, scaling='density')
f = (np.roll(f, len(f)/2))
pxx = np.roll(pxx, len(pxx)/2)
pxx = 10*np.log10(pxx)
plt.plot(f, pxx, label='Decimated')

# plot the low-pass filtered signal
f, pxx = scisig.welch(filtered_samples, fs=sample_rate, nperseg=nperseg, scaling='density')
f = (np.roll(f, len(f)/2))
pxx = np.roll(pxx, len(pxx)/2)
pxx = 10*np.log10(pxx)

plt.plot(f, pxx, label='original signal')
plt.xlim([-45e3, 45e3])
plt.ylim([np.min(full_pxx)-40, np.max(pxx)+1])

plt.title("Spectrum of one channel at baseband")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
plt.legend(loc='best')

# plot the signal raised to the 4th power
# Gives an idea of frequency offset after doppler compensation
plt.subplot(223)
f, pxx = scisig.welch(signal_to_4th_power, fs=sample_rate/decimation, nperseg=nperseg, scaling='density')
f = (np.roll(f, len(f)/2))
pxx = np.roll(pxx, len(pxx)/2)
pxx = 10*np.log10(pxx)

plt.plot(f, pxx)
plt.xlim([-2e3, 2e3])
plt.title("Spectrum of signal to 4th power")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")


# Plot complex samples from one channel
ax = plt.subplot(224)
start = int(1e6)
stop = start + int(100e3)
decim = 8*32
filtered_samples /= np.max(np.abs(filtered_samples))
baseband_samples /= np.max(np.abs(baseband_samples))
plt.scatter(filtered_samples.real[start:stop:decim], filtered_samples.imag[start:stop:decim])
plt.scatter(baseband_samples.real[int(start/decimation):int(stop/decimation):8], baseband_samples.imag[int(start/decimation):int(stop/decimation):8])
plt.title("Complex samples (10k samples)")
plt.xlabel("Real")
plt.ylabel("Imag")
ax.set_aspect(aspect=1)

plt.tight_layout()
plt.show()
