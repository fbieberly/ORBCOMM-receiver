##############################################################################
#
# Author: Frank Bieberly
# Date: 6 February 2020
# Name: record_orbcomm_long.py
# Description:
# This script is used to make longer recordings of an orbcomm signal. To reduce
# the size of the recordings, the signal is filtered and decimated down to
# 19.2 kHz.
#
#
##############################################################################

from CONFIG import lat, lon, alt, min_elevation

import signal
from time import time, sleep
from math import degrees, log10

# Need this to enable ctrl+c to exit the program.
signal.signal(signal.SIGINT, signal.SIG_DFL)

import ephem
import numpy as np
from rtlsdr import RtlSdr
from scipy.io import savemat
from scipy.signal import firwin

from helpers import get_tle_lines
from sat_db import active_orbcomm_satellites

# Create a pyephem sat object for all the active satellites
# using latest TLE data
for name in active_orbcomm_satellites:
    sat_line0, sat_line1, sat_line2 = get_tle_lines(name, tle_dir='./tles')
    sat = ephem.readtle(str(sat_line0), str(sat_line1), str(sat_line2))
    active_orbcomm_satellites[name]['sat_obj'] = sat
    active_orbcomm_satellites[name]['tles'] = [sat_line0, sat_line1, sat_line2]


# PyEphem observer
# lat/lon/alt defined in CONFIG
obs = ephem.Observer()
obs.lat, obs.lon = '{}'.format(lat), '{}'.format(lon)
obs.elevation = alt

# speed of light
c = 299792458.0 # m/s

# Set RTLSDR parameters and initialize
# This sample rate is used because it is a multiple of the baud rate
# Also it allows capturing the whole 1 MHz channel, which ensures getting
# both orbcomm channels without tuning to a different center frequency
# However, it is a much higher sample rate than would be needed.
sample_rate = 1.2288e6
decimation =  64.0      # brings sample rate to 19200.0 Hz
center_freq = 137.5e6   # we will change this later.
gain = 'auto' # Use AGC

sdr = RtlSdr()
sdr.rs = sample_rate
sdr.gain = gain
sdr.fc = center_freq


# receive samples that are an integer multiple of 1024 from the RTLSDR
num_samples_per_recording = int(1024*128)
should_finish = False
queue_max_size = 30
list_of_arrays = [] # Is this a horrible idea?
record_time = None
lpf_remaining_samples = np.array([])


# How long to record by
record_duration = 60.0     # seconds
needed_arrays = record_duration * sample_rate / num_samples_per_recording
max_arrays = int(needed_arrays)

# This is a callback function for async rtlsdr receive samples
def rtlsdr_callback(samples, context):
    global should_finish
    global list_of_arrays
    global lpf_remaining_samples
    global record_time
    global lat, lon, alt
    global max_arrays
    if record_time == None:
        record_time = time()

    sdr = context['sdr']
    sat = context['sat']
    sat_name = context['sat_name']
    obs = context['observer']
    tles = context['tles']
    center_freq = context['fc']

    # Initialize LPF parameters
    lpf_filter_freq = 10e3
    lpf_order = 181   # filter order
    lpf_nyq = 0.5 * sample_rate
    lpf_normal_cutoff = lpf_filter_freq / lpf_nyq
    lpf_fir_taps = firwin(lpf_order, lpf_filter_freq, fs=sample_rate)

    # Calculate the sat position
    obs.date = ephem.now()
    sat.compute(obs)

    # Normalize samples
    samples /= np.median(np.abs(samples))

    # save samples to .mat file
    sample_buffer = np.concatenate([lpf_remaining_samples, samples])

    samps_to_filter = int(len(sample_buffer) - lpf_order)
    while samps_to_filter % decimation != 0:
        samps_to_filter -= 1

    # Low pass filter and decimate in one step
    decimated_samples = np.zeros(int(samps_to_filter / decimation), dtype=np.complex64)
    for yy in range(0, len(decimated_samples)):
        temp_samples = sample_buffer[int(yy*decimation):int(yy * decimation + lpf_order)]
        decimated_samples[yy] = np.dot(lpf_fir_taps, temp_samples)
    lpf_remaining_samples = sample_buffer[int((yy + 1) * decimation):]

    list_of_arrays.append(decimated_samples)
    print("Length of list_of_arrays: {}".format(len(list_of_arrays)))

    # If the satellite goes below the horizon, stops receiving samples
    if degrees(sat.alt) < min_elevation or len(list_of_arrays) > max_arrays:

        big_array = np.concatenate(list_of_arrays)
        filename = '{:.3f}'.format(record_time).replace('.', 'p') + '.mat'
        save_dict = {
                    'samples':big_array,
                    'timestamp':record_time,
                    'sats': sat_name,
                    'tles': tles,
                    'fs': sample_rate/decimation,
                    'fc': center_freq,
                    'lat':lat,
                    'lon':lon,
                    'alt':alt,
                    }
        savemat('./data/' + filename, save_dict, do_compression=True)
        print("File saved: {}".format('./data/' + filename))
        should_finish = True
        sdr.cancel_read_async()
        return 0

# This is the main script loop
while 1:
    try:
        start_loop = time()
        obs.date = ephem.now()
        if should_finish:
            break

        # Find satellites above the horizon
        sats = []
        for sat_name in active_orbcomm_satellites:
            sat = active_orbcomm_satellites[sat_name]['sat_obj']
            sat.compute(obs)
            if degrees(sat.alt) > min_elevation:
                sats.append((sat_name, sat, degrees(sat.alt)))

        if len(sats) > 0:
            # Find the satellite that has the highest elevation.
            print("\nSatellites overhead: ")
            sorted_sats = sorted(sats, key=lambda x: x[2], reverse=True)
            for sat_name, sat, degrees_above_horizon in sorted_sats:
                print('{:20}: {:3.1f} degrees elevation'.format(sat_name, degrees_above_horizon))

            sat_name = sorted_sats[0][0]
            sat = sorted_sats[0][1]
            tles = active_orbcomm_satellites[sat_name]['tles']

            print("Receiving from: {}".format(sat_name))
            frequencies = active_orbcomm_satellites[sorted_sats[0][0]]['frequencies']
            print("Satellite frequencies: {}".format(frequencies))
            # Decode the lower of the two channels
            sat_center_frequency = frequencies[0]
            center_freq = sat_center_frequency   # we will change this later.
            sdr.fc = center_freq


            print('Recording samples.')
            # Record samples twice just to fill up buffers (not sure if needed)
            samples = sdr.read_samples(num_samples_per_recording)
            samples = sdr.read_samples(num_samples_per_recording)

            context_dict = {
                            'sdr':sdr,
                            'observer':obs,
                            'sat':sat,
                            'sat_name':sat_name,
                            'tles':[tles],
                            'fs': sample_rate/decimation,
                            'fc': center_freq,
                            }
            sdr.read_samples_async(rtlsdr_callback, num_samples_per_recording, context_dict)

            if should_finish:
                break

        else:
            # If no satellite is overhead, find the next one that will be
            sat_detected = False
            for minute in range(0, 60*12):

                obs.date = ephem.now() + minute * ephem.minute
                for sat_name in active_orbcomm_satellites:
                    sat = active_orbcomm_satellites[sat_name]['sat_obj']
                    sat.compute(obs)
                    if degrees(sat.alt) > min_elevation:
                        sat_detected = True

                        if minute > 1:
                            print("Time until next satellite  ({}) visible: {:.0f} minutes".format(sat_name, minute))
                            sleep(60)
                        else:
                            sleep(1)

                        break
                if sat_detected:
                    break
            if sat_detected == False:
                print("No upcoming satellite passes detected within 12 hours. Exiting.")
                exit()

    except KeyboardInterrupt:
        break

sdr.close()
print("Program exited.")