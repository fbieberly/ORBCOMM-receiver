##############################################################################
#
# Author: Frank Bieberly
# Date: 30 April 2019
# Name: file_decoder.py
# Description: 
# This script takes in .mat files (produced by record_spectrum_orbcomm.py) and
# produces multiple plots of the signals spectrum, constellation, timing 
# recovery, eye diagram, and IQ samples. Additionally, decoded bits are printed
# and saved to a file.
#
##############################################################################

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
from helpers import butter_lowpass_filter, complex_mix, rrcosfilter


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

# Normalize samples
samples /= np.median(np.abs(samples))

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
samples_per_symbol = 2  # We should only need 2 samples per symbol
decimation = int(sample_rate/(samples_per_symbol*baud_rate))
decimated_samples = filtered_samples[::decimation]


# estimate remaining carrier error (RTLSDR frequency error)
# signal to the fourth power, take the fft, peak is at frequency offset
rbw = 1
nperseg = int(sample_rate/decimation/rbw)
signal_to_4th_power = np.power(decimated_samples, 4)
f, pxx = scisig.welch(signal_to_4th_power, fs=sample_rate/decimation, \
                      return_onesided=False, nperseg=nperseg, scaling='density')
f = (np.roll(f, len(f)/2))
pxx = np.roll(pxx, len(pxx)/2)
search_window = int(1000. / ((sample_rate/decimation)/nperseg))  # search +/- 1 kHz around fc
frequency_peak = np.argmax(pxx[len(pxx)/2 - search_window:len(pxx)/2 + search_window])
freq_offset = -(frequency_peak - search_window)*(sample_rate/decimation/nperseg) / 4
baseband_samples = complex_mix(decimated_samples, freq_offset, sample_rate/decimation)
print "Remaining frequency offset: {}".format(freq_offset)

# Create RRC taps
alpha = 0.4
baudrate = 1.
num_of_symbols_half_filter = 8.
rrc_num_taps = samples_per_symbol * num_of_symbols_half_filter * 2. + 1.
t_idx, rrc_taps = rrcosfilter(rrc_num_taps, alpha, baudrate, samples_per_symbol)
matched_filtered_samples = scisig.lfilter(rrc_taps, [1.0], baseband_samples)

# Manually find a close timing offset
# For only 2 samples per symbol, 0 works all the time.
sample_delay = 0

tau     = 0.                        # initial timing offset estimate
dtau    = 0.                        # initial timing _rate_ offset estimate
buf     = np.zeros(5, dtype=np.complex64) # filter buffer
th      = np.array([-2., -1., 0., 1., 2.]) # time vector for interpolating filter
w       = np.sqrt(np.hamming(5).T)  # window function for interpolating filter
alpha   = 0.0005                    # loop filter bandwidth (timing _rate_ adjustment factor)
beta    = 2*np.sqrt(alpha)          # (timing _phase_ adjustment factor)
counter = sample_delay + 1          # interpolating filter adds delay
time_recovery_samples       = np.zeros(len(matched_filtered_samples), dtype=np.complex64)
buf_dz  = [0.,0.,0.]

dtau_vect = np.zeros(len(matched_filtered_samples), dtype=np.complex64)
tau_vect = np.zeros(len(matched_filtered_samples), dtype=np.complex64)

for i in range(1, len(matched_filtered_samples)):
    # push sample into interpolating filter
    buf[:-1] = buf[1:]
    buf[-1] = matched_filtered_samples[i]

    # interpolate matched filter output
    hi   = np.sinc(th - tau) * w  # interpolating filter coefficients
    hi /= np.sum(hi)
    time_recovery_samples[i] = np.dot(buf, np.flip(hi))    # compute matched filter output

    # take (approximate) derivative of filter output
    buf_dz[:-1] = buf_dz[1:]
    buf_dz[-1] = time_recovery_samples[i]
    dz = -np.dot(buf_dz, np.array([-1, 0, 1]))

    # determine if an output sample needs to be computed
    counter = counter + 1
    if counter >= samples_per_symbol:
        # decrement counter by samples per symbol
        counter = counter - samples_per_symbol

        # compute timing error signal, accounting for delay
        err = np.tanh( (dz * np.conj(time_recovery_samples[i-1])).real )

        # update timing estimate
        tau = tau + alpha*err

        # update timing rate change
        dtau = dtau + alpha * err
        tau  = tau  +  beta * err

    # update timing error
    tau = tau + dtau/samples_per_symbol

    # save results for plotting
    dtau_vect[i] = dtau
    tau_vect[i] = tau

# Plot timing offset
plt.figure()
plt.subplot(211)
plt.title("Timing offset (tau)")
plt.plot(tau_vect)

plt.subplot(212)
plt.title("Derivative (Dtau)")
plt.plot(dtau_vect)
plt.tight_layout()

# Plot eye-diagram
plt.figure()
plt.subplot(311)
plt.title("Eye Diagram before matched filter")
offset = samples_per_symbol * 200
num_plots = 8
length = 64
for xx in xrange(num_plots):
    plt.plot(baseband_samples[offset:offset+length].imag)
    offset += length
plt.subplot(312)
plt.title("After matched filter")

offset = samples_per_symbol * 200
for xx in xrange(num_plots):
    plt.plot(matched_filtered_samples[offset:offset+length].imag)
    offset += length
plt.grid()
plt.subplot(313)
plt.title("After timing recovery")

offset = samples_per_symbol * 200
for xx in xrange(num_plots):
    plt.plot(time_recovery_samples[offset:offset+length].imag)
    offset += length
plt.grid()
plt.tight_layout()

# After timing recovery, performing a fine-frequency PLL
# First take the signal to the fourth power, then low pass filter
filter_freq = 0.1e3
signal_to_4th_power2 = np.power(time_recovery_samples, 4)
filtered_sig4th = butter_lowpass_filter(signal_to_4th_power2, filter_freq, sample_rate/decimation, order=5)

# PLL code loosely based on liquid dsp simple pll tutorial: 
# http://liquidsdr.org/blog/pll-simple-howto/
phase_est = np.zeros(len(time_recovery_samples)+1)
alpha = 0.05
beta = 0.5 * alpha**2
frequency_out = 0.0
phase_out = []

for idx, sample in enumerate(filtered_sig4th):
    signal_out = np.exp(1j * phase_est[idx])
    phase_error = np.angle( sample * np.conj(signal_out) )
    phase_est[idx+1] = phase_est[idx] + alpha * phase_error
    frequency_out += beta * phase_error
    phase_out.append(frequency_out)


# Phase compensate the IQ samples
phase_comp_samples = time_recovery_samples * np.conj(np.exp(1j*phase_est[:-1]/4.)) * np.conj(np.exp(1j*np.pi/4.))

# Decode to bits
# normalize phase compensated samples;
phase_comp_samples /= np.median(np.abs(phase_comp_samples))

# Select peak sample from each symbol
demod_symbols = phase_comp_samples[::samples_per_symbol]

#Differential demodulation
# A 1 is a 90 degree phase shift forward from the last symbol
# A 0 is a -90 degree phase shift forward from the last symbol
bits = []
angles = []
for xx in range(1, len(demod_symbols)):
    angle = np.angle(demod_symbols[xx], deg=True) - np.angle(demod_symbols[xx-1], deg=True)
    if angle > 180:
        angle -= 360
    if angle < -180:
        angle += 360
    angles.append(angle)
    bit = 0
    if angle > 0: bit = 1
    bits.append(bit)

# Plot the angles between each symbol, should be +/- 90 degrees
plt.figure()
plt.title('Angle between successive symbols')
plt.xlabel('symbol number')
plt.ylabel('Angle (degrees)')
plt.plot(angles, 'x')
plt.grid()

# The bits are xor'ed, so xor them again to retreive the original bit stream
xord_bits = []
for xx in range(1, len(bits)):
    xor = 0
    if np.abs(sum(bits[xx-1:xx])) == 1.0:
        xor = 1
    xord_bits.append(xor)

#print out the bits
print("Bits:")
print(''.join([str(bit) for bit in xord_bits]))


# Plot IQ samples
plt.figure()
plt.subplot(211)
plt.title("After carrier recovery")
plt.plot(phase_comp_samples[1000:1080].real)
plt.plot(phase_comp_samples[1000:1080].imag)
plt.grid()

plt.subplot(212)
plt.title("After carrier recovery")
plt.plot(time_recovery_samples[1000:1080].real)
plt.plot(time_recovery_samples[1000:1080].imag)
plt.grid()
plt.tight_layout()


# Plot spectrum of recording
plt.figure()
plt.subplot(221)
nperseg = int(sample_rate/100.0)
f, full_pxx = scisig.welch(samples, fs=sample_rate, nperseg=nperseg, \
                           return_onesided=False, scaling='density')
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
f, pxx = scisig.welch(baseband_samples, fs=sample_rate/decimation, \
                      return_onesided=False, nperseg=nperseg, scaling='density')
f = (np.roll(f, len(f)/2))
pxx = np.roll(pxx, len(pxx)/2)
pxx = 10*np.log10(pxx)
plt.plot(f, pxx, label='Decimated')

# plot the low-pass filtered signal
f, pxx = scisig.welch(filtered_samples, fs=sample_rate, nperseg=nperseg, \
                      return_onesided=False, scaling='density')
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
nperseg = len(signal_to_4th_power)
f, pxx = scisig.welch(signal_to_4th_power, fs=sample_rate/decimation, \
                       return_onesided=False, nperseg=nperseg, scaling='density')
f = (np.roll(f, len(f)/2))
pxx = np.roll(pxx, len(pxx)/2)
pxx = 10*np.log10(pxx)

plt.plot(f, pxx)
plt.xlim([-2e3, 2e3])
plt.title("Spectrum of signal to 4th power")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
plt.tight_layout()


# Plot complex samples from one channel
plt.figure()
ax = plt.subplot(111)

matched_filtered_samples /= np.median(np.abs(matched_filtered_samples))
phase_comp_samples /= np.median(np.abs(phase_comp_samples))

plt.scatter(matched_filtered_samples.real[100*samples_per_symbol::samples_per_symbol], matched_filtered_samples.imag[100*samples_per_symbol::samples_per_symbol], marker='x', label='after MF')
plt.scatter(phase_comp_samples.real[100*samples_per_symbol::samples_per_symbol], phase_comp_samples.imag[100*samples_per_symbol::samples_per_symbol], marker='x', label='timing recovery')
plt.legend(loc='best')
plt.title("Complex samples (10k samples)")
plt.xlabel("Real")
plt.ylabel("Imag")
ax.set_aspect(aspect=1)
plt.tight_layout()
plt.show()
