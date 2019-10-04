##############################################################################
#
# Author: Frank Bieberly
# Date: 17 July 2019
# Name: realtime_receiver_tkgui.py
# Description:
# This script uses the realtime decoder class to decode one data channel from
# an orbcomm satellite that is overhead. When the satellite goes below the
# horizon, the script will search for the next satellite with the highest
# elevation. If none is overhead, it will wait until one is.
# The GUI has a few slightly useful bits of information for troubleshooting
# if you aren't getting a good reception.
# The code uses multiprocessing to run the decoder in a separate thread from
# the RTLSDR code. I haven't gotten around to getting the code to finish
# cleanly.
#
#
##############################################################################


import matplotlib
qt_backend_found = False

for backend in ['Qt4Agg', 'Qt5Agg']:
    try:
        matplotlib.use(backend)
        qt_backend_found = True
        break
    except ValueError:
        pass
if qt_backend_found == False:
    print("Qt plotting backend for matplotlib required. Exiting.")
    exit()

import signal
import cPickle as pickle
from time import time, sleep
import multiprocessing as mp
from math import degrees, log10

# Need this to enable ctrl+c to exit the program.
signal.signal(signal.SIGINT, signal.SIG_DFL)

import ephem
import numpy as np
from rtlsdr import RtlSdr
from scipy.signal import welch
import matplotlib.pyplot as plt


from helpers import get_tle_lines
from realtime_decoder_new import RealtimeDecoder
from sat_db import active_orbcomm_satellites

from CONFIG import *

# Create a pyephem sat object for all the active satellites
# using latest TLE data
for name in active_orbcomm_satellites:
    sat_line0, sat_line1, sat_line2 = get_tle_lines(name, tle_dir='./tles')
    sat = ephem.readtle(sat_line0, sat_line1, sat_line2)
    active_orbcomm_satellites[name]['sat_obj'] = sat

# PyEphem observer
# lat/lon/alt defined in CONFIG
obs = ephem.Observer()
obs.lat, obs.lon = '{}'.format(lat), '{}'.format(lon)
obs.elevation = alt

##########################################
# configure the matplotlib plots
img = plt.imread("./map.jpg")
const_fig, ax_arr = plt.subplots(2, 2, num=10)

# Plot of the complex IQ samples
cost_line, = ax_arr[0,0].plot([0,], 'bx')
ax_arr[0,0].set_xlim((-2, 2))
ax_arr[0,0].set_ylim((-2, 2))
ax_arr[0,0].set_title("Complex I/Q samples")
ax_arr[0,0].set_xlabel("I")
ax_arr[0,0].set_ylabel("Q")

# Plot of the channel spectrum
fft_line, = ax_arr[1,0].plot([0,], 'b')
ax_arr[1,0].set_xlim((-6e3, 6e3))
ax_arr[1,0].set_ylim((-50, 20))
ax_arr[1,0].set_title("Channel Spectrum")
ax_arr[1,0].set_xlabel("Frequency (centered on channel) Hz")
ax_arr[1,0].set_ylabel("Power (dB)")

# Plot of globe and ground station location
img_plot = ax_arr[1,1].imshow(img, extent=[-180, 180, -90, 90])
loc_plot, = ax_arr[1,1].plot(lon, lat, 'ro')

ax_arr[0,1].axis('off')

##########################################

# speed of light
c = 299792458.0 # m/s

# Set RTLSDR parameters and initialize
# This sample rate is used because it is a multiple of the baud rate
# Also it allows capturing the whole 1 MHz channel, which ensures getting
# both orbcomm channels without tuning to a different center frequency
# However, it is a much higher sample rate than would be needed.
sample_rate = 1.2288e6
center_freq = 137.5e6
gain = 'auto' # Use AGC

sdr = RtlSdr()
sdr.rs = sample_rate
sdr.gain = gain
sdr.fc = center_freq

# receive samples that are an integer multiple of 1024 from the RTLSDR
num_samples_per_recording = int(1024*128)
should_finish = False

def rtlsdr_callback(samples, context):
    global decoder
    global should_finish

    sdr = context['sdr']
    sat = context['sat']
    sat_name = context['sat_name']
    obs = context['observer']

    obs.date = ephem.now()
    sat.compute(obs)
    relative_vel = sat.range_velocity
    doppler = c/(c+relative_vel) * sat_center_frequency - sat_center_frequency

    # This code will catch the case when the user closes the matplotlib figure
    if queue.qsize() > 20:
        should_finish = True
        queue.put((None, None))
        sdr.cancel_read_async()
        return 0

    # If the satellite goes below the horizon, stops receiving samples
    if degrees(sat.alt) < min_elevation:
        queue.put((None, None))
        sdr.cancel_read_async()
        return 0

    # I pass a new contex to the processing thread
    new_context = {
            'doppler':doppler,
            'sat_name':sat_name,
            }
    queue.put((samples, new_context))


def process_samples(queue):
    global decoder
    global should_finish
    tic = time()

    sat_gps_dict = {}
    sat_plot_lines = []

    while 1:
        try:
            samples, context = queue.get()
            if context == None:
                break

            doppler = context['doppler']
            sat_name = context['sat_name']

            tic = time()
            packet_list = decoder.decode_samples_to_packets(samples, doppler)
            # if should_save_to_file:
            #     do saving stuff
            decoder.parse_packets()

            # Get the welch spectrum for just the channel we are decoding
            f, pxx = welch(decoder.decimated_samples, fs=decoder.sample_rate/decoder.decimation, \
                           return_onesided=False, nperseg= int(len(decoder.decimated_samples)/10), \
                           scaling='density')
            f = (np.roll(f, int(len(f)/2)))
            pxx = np.roll(pxx, int(len(pxx)/2))
            pxx = 20*np.log10(pxx)
            pxx -= np.median(pxx)

            symbols = decoder.symbols
            symbols /= np.median(np.abs(symbols))

            # plot IQ samples
            cost_line.set_xdata(symbols.real)
            cost_line.set_ydata(symbols.imag)
            ax_arr[0,0].draw_artist(ax_arr[0,0].patch)
            ax_arr[0,0].draw_artist(cost_line)

            # plot spectrum
            fft_line.set_xdata(f)
            fft_line.set_ydata(pxx)
            ax_arr[1,0].draw_artist(ax_arr[1,0].patch)
            ax_arr[1,0].draw_artist(fft_line)

            # Check if we have new lat/long for satellite
            # If we do, plot the track of the satellite since we've been recording
            # This track is from received lat/long (not from TLE position)
            if decoder.sat_lon != 0.0 and decoder.sat_lat != 0.0:
                if abs(decoder.sat_lon) < 180 and abs(decoder.sat_lat) < 90:
                    sat_lat_lon = (decoder.sat_lon, decoder.sat_lat)
                    if sat_name not in sat_gps_dict:
                        sat_gps_dict[sat_name] = [sat_lat_lon]
                        sat_plot_lines.append(ax_arr[1,1].plot([0,], 'D-', color='cyan', markevery=[-1])[0])
                    else:
                        if sat_lat_lon not in sat_gps_dict[sat_name]:
                            sat_gps_dict[sat_name].append(sat_lat_lon)

                    # Only update the plot when we get a new lat/long
                    ax_arr[1,1].draw_artist(ax_arr[1,1].patch)
                    ax_arr[1,1].draw_artist(img_plot)
                    ax_arr[1,1].draw_artist(loc_plot)
                    for idx, key in enumerate(sat_gps_dict):
                        lon_lat_arr = sat_gps_dict[key]
                        lons, lats = zip(*lon_lat_arr)
                        sat_plot_lines[idx].set_xdata(lons)
                        sat_plot_lines[idx].set_ydata(lats)
                        text = ax_arr[1,1].text(lons[-1], lats[-1], "{}".format(key), color='cyan', wrap=True)
                        ax_arr[1,1].draw_artist(sat_plot_lines[idx])
                        ax_arr[1,1].draw_artist(text)

            toc = time()
            # The time ratio is the amount of time it takes to process a batch of samples
            # divided by the time that it takes to record the samples.
            # A time ratio over 1.0 means that you are falling behind.
            time_ratio = sdr.sample_rate/len(samples) * (toc-tic)

            bad_percent = decoder.bad_packets / (decoder.bad_packets + decoder.good_packets + 0.1)
            good_percent = decoder.good_packets / (decoder.bad_packets + decoder.good_packets + 0.1)

            text0 = ax_arr[0,1].text(0.10, 0.90, "Satellite Name: {}".format(sat_name))
            text1 = ax_arr[0,1].text(0.10, 0.80, "Bad Packets %:  {:4.1f}".format(bad_percent*100.0))
            text2 = ax_arr[0,1].text(0.10, 0.70, "Good Packets %: {:4.1f}".format(good_percent*100.0))
            text3 = ax_arr[0,1].text(0.10, 0.60, "Time ratio:     {:4.2f}".format(time_ratio))

            ax_arr[0,1].draw_artist(ax_arr[0,1].patch)
            ax_arr[0,1].draw_artist(text0)
            ax_arr[0,1].draw_artist(text1)
            ax_arr[0,1].draw_artist(text2)
            ax_arr[0,1].draw_artist(text3)

            const_fig.canvas.update()
            const_fig.canvas.flush_events()

            # If the user closes the figure, fignums goes to 0
            # Then we'll close this process and exit the program
            fignums = plt.get_fignums()
            if len(fignums) == 0:
                should_finish = True
                break

        # Catch ctrl+c and exit program
        except KeyboardInterrupt:
            break
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
            frequencies = active_orbcomm_satellites[sorted_sats[0][0]]['frequencies']
            print("Satellite frequencies: {}".format(frequencies))
            # Decode the lower of the two channels
            sat_center_frequency = frequencies[0]
            decoder = RealtimeDecoder(sat_center_frequency, center_freq=center_freq, sample_rate=sample_rate)

            # Calculate the starting doppler
            obs.date = ephem.now()
            sat.compute(obs)
            relative_vel = sat.range_velocity
            doppler = c/(c+relative_vel) * sat_center_frequency - sat_center_frequency

            # Open matplotlib figure for plotting:
            plt.pause(0.0001)

            print('Recording samples.')
            # Record samples twice just to fill up buffers (not sure if needed)
            samples = sdr.read_samples(num_samples_per_recording)
            samples = sdr.read_samples(num_samples_per_recording)

            print("Create multiprocessing queue")
            queue = mp.Queue()

            print("Creating multiprocessing process")
            p = mp.Process(target=process_samples, args=(queue,))
            p.start()

            # This code helps the decoder determine frequency offsets, etc.
            decoder.first_samples(samples, doppler=doppler)

            context_dict = {
                            'sdr':sdr,
                            'observer':obs,
                            'sat':sat,
                            'sat_name':sat_name,
                            }
            sdr.read_samples_async(rtlsdr_callback, num_samples_per_recording, context_dict)
            print("Ending async RTLSDR processing")
            queue.close()
            queue.join_thread()
            p.join()

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