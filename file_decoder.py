##############################################################################
#
# Author: Frank Bieberly
# Date: 30 April 2019
# Name: file_decoder.py
# Description:
# This script takes in .mat files (produced by record_spectrum_orbcomm.py) and
# produces multiple plots of the signals spectrum, constellation, timing
# recovery, eye diagram, and IQ samples. Additionally, raw decoded bits are
# printed and saved to a file (./packets.txt).
#
##############################################################################

import glob
from datetime import datetime, timedelta

import ephem
import numpy as np
from scipy.io import loadmat
import scipy.signal as scisig
import matplotlib.pyplot as plt

from sat_db import active_orbcomm_satellites
from orbcomm_packet import packet_dict
from helpers import butter_lowpass_filter, complex_mix, rrcosfilter
from helpers import fletcher_checksum, ecef_to_lla


# speed of light
c = 299792458.0 # m/s

# Perform brute force error correction
brute_force_1bit_error_correction = True

# Where to save decoded packets
packet_file = r'./packets.txt'

# Where the data files are located
data_dir = r'./data/'
sample_file = sorted(glob.glob(data_dir + "*.mat"))[-4]

# Load the .mat file and print some of the metadata
data = loadmat(sample_file)
print("Filename: {}".format(sample_file))
print("Timestamp: {}".format(data['timestamp'][0][0]))
print("Data collected on: {}".format(datetime.utcfromtimestamp(data['timestamp'][0][0])))
print("Satellites in recording: {}".format(', '.join(data['sats'])))
print("SDR Sample rate: {} Hz".format(data['fs'][0][0]))
print("SDR Center frequency: {} Hz".format(data['fc'][0][0]))
frequencies = []
for sat_name in data['sats']:
    freq1, freq2 = active_orbcomm_satellites[sat_name]['frequencies']
    frequencies.append(freq1)
    frequencies.append(freq2)
print("Satellite frequencies: {}".format(', '.join([str(xx) for xx in frequencies])))

# Extract the values from file for some further processing
samples = data['samples'][0]
center_freq = data['fc'][0][0]
sample_rate = data['fs'][0][0]
timestamp = data['timestamp'][0][0]
lat = data['lat'][0][0]
lon = data['lon'][0][0]
alt = data['alt'][0][0]
print("Number of samples: {}".format(len(samples)))

# Check which satellite frequency is contained in the recording
# If both frequencies are present, decode the lower one.
sat_center_frequency = 0
for freq in frequencies:
    if freq > (center_freq - sample_rate/2.0) and freq < (center_freq + sample_rate/2.0):
        sat_center_frequency = freq
        break
if sat_center_frequency == 0:
    print("Satellite channel frequencies are not in the recorded spectrum.")
    exit()
# To force decoding the upper frequency uncomment this line
# sat_center_frequency = frequencies[1]


# PyEphem observer
obs = ephem.Observer()
obs.lat, obs.lon = '{}'.format(lat), '{}'.format(lon)
obs.elevation = alt # Technically is the altitude of observer
obs.date = datetime.utcfromtimestamp(timestamp)

# Normalize samples
samples /= np.median(np.abs(samples))

# Get the TLE information from the .mat file
sat_line0, sat_line1, sat_line2 = [str(xx) for xx in data['tles'][0]]
sat = ephem.readtle(sat_line0, sat_line1, sat_line2)
sat.compute(obs)

# Use the TLE info that was in the .mat file to calculate doppler shift
# of the satellite's transmission
relative_vel = sat.range_velocity
doppler = c/(c+relative_vel) * sat_center_frequency - sat_center_frequency

# Mix the samples to baseband (compensating for doppler shift)
# There will be a residual carrier error because of the RLTSDR frequency offset
freq_shift = center_freq - sat_center_frequency - doppler
mixed_down_samples, _ = complex_mix(samples,
                                    freq_shift,
                                    sample_rate)

filter_freq = 10e3
baud_rate = 4800.0
samples_per_symbol = 2  # We should only need 2 samples per symbol

if sample_rate == 1.2288e6:
    # Low pass filter the signal
    filtered_samples = butter_lowpass_filter(mixed_down_samples, filter_freq, sample_rate, order=5)

    # Decimated signal
    decimation = int(sample_rate/(samples_per_symbol*baud_rate))
    decimated_samples = filtered_samples[::decimation]

elif sample_rate == 19200.0:
    print("Single channel recording detected.")
    filter_freq = 4e3
    filtered_samples = butter_lowpass_filter(mixed_down_samples, filter_freq, sample_rate, order=2)

    decimation = 2
    decimated_samples = filtered_samples[::decimation]
    filtered_samples = filtered_samples


# estimate remaining carrier error (RTLSDR frequency error)
# signal to the fourth power, take the fft, peak is at frequency offset
rbw = 1
nperseg = int(sample_rate/decimation/rbw)

# take first min(50e3, len(samples)) of the samples to the 4th power.
signal_to_4th_power = np.power(decimated_samples[:int(min(len(decimated_samples), 50e3)):], 4)
f, pxx = scisig.welch(signal_to_4th_power,
                      fs=sample_rate/decimation,
                      return_onesided=False,
                      nperseg=nperseg,
                      scaling='density')
f = (np.roll(f, int(len(f)/2)))
pxx = np.roll(pxx, int(len(pxx)/2))
search_window = int(1000. / ((sample_rate/decimation)/nperseg))  # search +/- 1 kHz around fc
frequency_peak = np.argmax(pxx[int(len(pxx)/2 - search_window):int(len(pxx)/2 + search_window)])
freq_offset = -(frequency_peak - search_window)*(sample_rate/decimation/nperseg) / 4
baseband_samples, _ = complex_mix(decimated_samples, freq_offset, sample_rate/decimation)
print("Remaining frequency offset after doppler compensation: {} Hz".format(freq_offset))

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

dtau_vect = np.zeros(len(matched_filtered_samples))
tau_vect = np.zeros(len(matched_filtered_samples))
timing_error = np.zeros(len(matched_filtered_samples))
err = 0.0

# Timing recovery
i = 1
i_offset = 0
while i < len(time_recovery_samples) - 1 and i - i_offset + 2 < len(matched_filtered_samples):
    # push sample into interpolating filter
    # If the time offset exceeds one sample in either direction
    if tau > 1.0:
        # Use the samples from last time
        i_offset += 1
        tau += -1
    elif tau < -1.0:
        # Push two samples into the buffer
        buf[:-2] = buf[2:]
        buf[-2:] = matched_filtered_samples[int(i - i_offset):int(i - i_offset + 2)]
        i_offset += -1
        tau += 1
    else:
        # Or just normally add one sample to the buffer
        buf[:-1] = buf[1:]
        buf[-1] = matched_filtered_samples[int(i - i_offset)]

    # interpolate matched filter output
    hi   = np.sinc(th - tau) * w  # interpolating filter coefficients
    hi /= np.sum(hi)
    time_recovery_samples[i] = np.dot(buf, np.flip(hi, 0))    # compute matched filter output

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

        # update timing rate change
        dtau = dtau + alpha * err
        tau  = tau  +  beta * err

    timing_error[i] = err
    # update timing error
    tau = tau + dtau/samples_per_symbol

    # save results for plotting
    dtau_vect[i] = dtau
    tau_vect[i] = tau + i_offset
    i += 1

# Plot timing offset
plt.figure()
plt.subplot(311)
plt.title("Timing offset (tau)")
plt.plot(tau_vect)

plt.subplot(312)
plt.title("Derivative (Dtau)")
plt.plot(dtau_vect)
plt.tight_layout()

plt.subplot(313)
plt.title("Timing error signal")
plt.plot(timing_error)
plt.tight_layout()

# Plot eye-diagram
plt.figure()
plt.subplot(311)
plt.title("Eye Diagram before matched filter")
offset = samples_per_symbol * 200
num_plots = 8
length = 64
for xx in range(num_plots):
    plt.plot(baseband_samples[offset:offset+length].imag)
    offset += length
plt.subplot(312)
plt.title("After matched filter")

offset = samples_per_symbol * 200
for xx in range(num_plots):
    plt.plot(matched_filtered_samples[offset:offset+length].imag)
    offset += length
plt.grid()
plt.subplot(313)
plt.title("After timing recovery")

offset = samples_per_symbol * 200
for xx in range(num_plots):
    plt.plot(time_recovery_samples[offset:offset+length].imag)
    offset += length
plt.grid()
plt.tight_layout()

# Costas loop
phase_est = np.zeros(len(time_recovery_samples)+1)
alpha = 0.03
beta = 0.2 * alpha**2
frequency_out = 0.0
freq_out = []

for idx, sample in enumerate(time_recovery_samples):
    signal_out = sample * np.exp(-1j * phase_est[idx])
    phase_error = np.sign(signal_out.real)*signal_out.imag - np.sign(signal_out.imag)*signal_out.real
    frequency_out += beta * phase_error
    phase_est[idx+1] = phase_est[idx] + alpha * phase_error + frequency_out
    freq_out.append(frequency_out)

# Alternative phase correcting code:
# PLL code loosely based on liquid dsp simple pll tutorial:
# http://liquidsdr.org/blog/pll-simple-howto/
## After timing recovery, performing a fine-frequency PLL
## First take the signal to the fourth power, then low pass filter
# filter_freq = 0.5e3
# signal_to_4th_power2 = np.power(time_recovery_samples, 4)
# filtered_sig4th = butter_lowpass_filter(signal_to_4th_power2, filter_freq, sample_rate/decimation, order=5)
# phase_est = np.zeros(len(time_recovery_samples)+1)
# alpha = 0.01
# beta = 0.2 * alpha**2
# frequency_out = 0.0
# freq_out = []
# for idx, sample in enumerate(filtered_sig4th):
#     signal_out = np.exp(1j * phase_est[idx])
#     phase_error = np.angle( sample * np.conj(signal_out) )
#     old_freq = frequency_out
#     frequency_out += beta * phase_error
#     phase_est[idx+1] = phase_est[idx] + alpha * phase_error + frequency_out
#     freq_out.append(frequency_out)

plt.figure()
plt.subplot(211)
plt.title('Phase output of PLL')
plt.plot(phase_est)
plt.grid()

plt.subplot(212)
plt.title('Frequency of PLL')
plt.plot(freq_out)
plt.grid()

# Phase compensate the IQ samples
phase_comp_samples = time_recovery_samples * np.conj(np.exp(1j*phase_est[:-1]))

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

# convert to a string of 0's and 1's
bit_string = ''.join([str(bit) for bit in bits])

# Find first full packet
size_of_packets = 12*8 # bits
num_of_possible_packets = len(bit_string)/size_of_packets
bit_offset = 0
print("Number of possible packets: {}".format(num_of_possible_packets))

# Extract the headers from the packet dictionary
packet_headers = [packet_dict[packet_type]['header'] for packet_type in packet_dict]

# for all bit offsets (of the length of the packets)
# calculate a score for most valid headers of that offset
# this also checks a bit-reversed score (in case my endianness is reversed)
scores = np.zeros(size_of_packets)
revscores = np.zeros(size_of_packets)
for xx in range(0, size_of_packets):
    for yy in range(xx, len(bit_string)-xx-8, size_of_packets):
        if bit_string[yy:yy+8][::-1] in packet_headers:
            scores[xx] += 1
        if bit_string[yy:yy+8] in packet_headers:
            revscores[xx] += 1

reverse = False
if np.max(scores) < np.max(revscores):
    reverse = True
if reverse:
    bit_offset = np.argmax(revscores)
else:
    bit_offset = np.argmax(scores)
print("Bit stream offset: {}".format(bit_offset))


packets = []
last_packet_epheris = False
for xx in range(bit_offset, len(bit_string)-size_of_packets, size_of_packets):
    if last_packet_epheris == True:
        last_packet_epheris = False
        continue
    packet = ''
    if reverse:
        header = '{:02X}'.format(int(bit_string[xx:xx+8], 2))
    else:
        header = '{:02X}'.format(int(bit_string[xx:xx+8][::-1], 2))

    ephemeris_header = packet_dict['Ephemeris']['hex_header']
    packet_length = 12*8
    if header == ephemeris_header:
        packet_length = 24*8
        last_packet_epheris = True

    for yy in range(0, packet_length, 8):
        if reverse:
            packet += '{:02X}'.format(int(bit_string[xx+yy:xx+yy+8], 2))
        else:
            packet += '{:02X}'.format(int(bit_string[xx+yy:xx+yy+8][::-1], 2))
    packets.append(packet)

# Save the packets (in hex) to a file
with open(packet_file, 'w') as f:
    for packet in packets:
        f.write(packet + '\n')

# Print out the parsed packets
error_packets = 0
fixed_packets = 0
print("\nList of packets: (### indicates checksum failed)")
for packet in packets:
    output = ''

    # Compute the fletcher16 checksum over the whole packet
    # 0000 output is a good packet
    if fletcher_checksum(packet) != '0000' and brute_force_1bit_error_correction:
        # Attempt to correct errors by flipping bits until the checksum comes out correct
        binary_packet = format(int(packet, 16), '0>{}b'.format(len(packet)/2*8))
        for xx in range(0, len(binary_packet)):
            temp_packet = ''
            if binary_packet[xx] == '0':
                temp_bits = binary_packet[:xx] + '1' + binary_packet[xx+1:]
            else:
                temp_bits = binary_packet[:xx] + '0' + binary_packet[xx+1:]

            for yy in range(0, len(binary_packet), 8):
                if reverse:
                    temp_packet += '{:02X}'.format(int(temp_bits[yy:yy+8][::-1], 2))
                else:
                    temp_packet += '{:02X}'.format(int(temp_bits[yy:yy+8], 2))
            if fletcher_checksum(temp_packet) == '0000':
                # print("Found correct packet!")
                packet = temp_packet
                fixed_packets += 1
                break

    if fletcher_checksum(packet) != '0000':
        output += '### '
        error_packets += 1
    for packet_type in packet_dict:
        packet_info = packet_dict[packet_type]
        if packet[:2] == packet_info['hex_header']:
            output += '{}: '.format(packet_type)
            for part, (start, stop) in packet_info['message_parts']:
                output += '{}: {} '.format(part, packet[start:stop])
            print(output)
            if packet_type == 'Ephemeris':
                payload = ''.join([packet[xx:xx+2] for xx in range(42, 2, -2)])

                # calculate current satellite time
                start_date = datetime.strptime('Jan 6 1980 00:00', '%b %d %Y %H:%M')
                week_number = payload[:4]
                time_of_week = payload[4:10]
                this_week = start_date + timedelta(weeks=int(week_number, 16))
                this_time = this_week + timedelta(seconds=int(time_of_week, 16))
                print("\tCurrent satellite time: {} Z".format(this_time))

                # calculate satellite ECEF position
                zdot = payload[10:15][::-1]
                ydot = payload[15:20][::-1]
                xdot = payload[20:25][::-1]
                zpos = payload[25:30][::-1]
                ypos = payload[30:35][::-1]
                xpos = payload[35:40][::-1]

                max_r_sat = 8378155.0
                val_20_bits = 1048576.0

                x_temp = int(xpos[:2][::-1], 16) + 256. * int(xpos[2:4][::-1], 16) + 256**2 * int(xpos[4:], 16)
                x_ecef = ((2*x_temp*max_r_sat)/val_20_bits - max_r_sat)
                y_temp = int(ypos[:2][::-1], 16) + 256. * int(ypos[2:4][::-1], 16) + 256**2 * int(ypos[4:], 16)
                y_ecef = ((2*y_temp*max_r_sat)/val_20_bits - max_r_sat)
                z_temp = int(zpos[:2][::-1], 16) + 256. * int(zpos[2:4][::-1], 16) + 256**2 * int(zpos[4:], 16)
                z_ecef = ((2*z_temp*max_r_sat)/val_20_bits - max_r_sat)

                lat, lon, alt = ecef_to_lla(x_ecef, y_ecef, z_ecef)
                print("\tLat/Lon:       {:8.4f}, {:8.4f}, Altitude: {:6.1f} km".format(lat, lon, alt/1000.0))
                print("\tEphem Lat/Lon: {:8.4f}, {:8.4f}, Altitude: {:6.1f} km".format(np.degrees(sat.sublat), np.degrees(sat.sublong), sat.elevation/1000.0))
            break
    # Unrecognized just means I don't know what these packets are for
    # would also happen if the header is corrupted
    if output in ['', '### ']:
        print("{}Unrecognized packet: {}".format(output, packet))

print("{} packets with errors, {} packets corrected, PER: {:5.1f}%".format(error_packets+fixed_packets,
                                                                           fixed_packets,
                                                                           float(error_packets)/len(packets)*100))

# Plot IQ samples
plt.figure()
plt.subplot(211)
plt.title("Before carrier recovery")
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
f = (np.roll(f, int(len(f)/2)) + center_freq)/1e6
full_pxx = np.roll(full_pxx, int(len(full_pxx)/2))
full_pxx = 10*np.log10(full_pxx)
plt.plot(f, full_pxx)
plt.title("Periodogram of recording\nRed dots are orbcomm channels")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
low_point = np.min(full_pxx)
for freq in frequencies:
    plt.plot(freq/1e6, low_point, 'ro')

# Plot one of the channels as baseband
# and plot the decimated signal
plt.subplot(222)
f, pxx = scisig.welch(baseband_samples, fs=sample_rate/decimation, \
                      return_onesided=False, nperseg=nperseg, scaling='density')
f = (np.roll(f, int(len(f)/2)))
pxx = np.roll(pxx, int(len(pxx)/2))
pxx = 10*np.log10(pxx)
plt.plot(f, pxx, label='Decimated')

# plot the low-pass filtered signal
f, pxx = scisig.welch(filtered_samples, fs=sample_rate, nperseg=nperseg, \
                      return_onesided=False, scaling='density')
f = (np.roll(f, int(len(f)/2)))
pxx = np.roll(pxx, int(len(pxx)/2))
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
f = (np.roll(f, int(len(f)/2)))
pxx = np.roll(pxx, int(len(pxx)/2))
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
plt.grid()
ax.set_aspect(aspect=1)
plt.tight_layout()
plt.show()
