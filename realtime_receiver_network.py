##############################################################################
#
# Author: Frank Bieberly
# Date: 17 July 2019
# Name: realtime_receiver_network.py
# Description:
# This script uses the realtime decoder class to decode one data channel from
# an orbcomm satellite that is overhead. When the satellite goes below the
# horizon, the script will search for the next satellite with the highest
# elevation. If none is overhead, it will wait until one is.
# The GUI has a few slightly useful bits of information for troubleshooting
# if you aren't getting a good reception.
# This script allows you to use a remote RTLSDR on a different computer (e.g.
# raspberry pi). The remote computer uses rtlsdr and netcat to send samples to
# the processing computer over port 5556.
# Usage:
# On processing computer: python realtime_receiver_network.py
# On remote rtlsdr computer: rtl_sdr -f 137500000 -s 1228800 -g 0 - | netcat -uv your_ip_address 5556
#
# your_ip_address is the IP address of the processing computer.
#
#
##############################################################################

from CONFIG import lat, lon, alt, min_elevation, realtime_plotting, ip_address
from CONFIG import port_number
if realtime_plotting:
    import matplotlib
    qt_backend_found = False

    for backend in ['Qt5Agg', 'Qt4Agg']:
        try:
            matplotlib.use(backend)
            qt_backend_found = True
            break
        except ValueError:
            pass
    if qt_backend_found == False:
        print("Qt plotting backend for matplotlib required. Exiting.")
        exit()

try:
    # For python2.7
    import cPickle as pickle
except:
    pass

import socket
import signal
from time import time, sleep
import multiprocessing as mp
from math import degrees, log10

# Need this to enable ctrl+c to exit the program.
signal.signal(signal.SIGINT, signal.SIG_DFL)

import ephem
import numpy as np
from scipy.signal import welch
import matplotlib.pyplot as plt

from helpers import get_tle_lines
from realtime_decoder import RealtimeDecoder
from sat_db import active_orbcomm_satellites

# Create a pyephem sat object for all the active satellites
# using latest TLE data
for name in active_orbcomm_satellites:
    sat_line0, sat_line1, sat_line2 = get_tle_lines(name, tle_dir='./tles')
    sat = ephem.readtle(str(sat_line0), str(sat_line1), str(sat_line2))
    active_orbcomm_satellites[name]['sat_obj'] = sat

# PyEphem observer
# lat/lon/alt defined in CONFIG
obs = ephem.Observer()
obs.lat, obs.lon = '{}'.format(lat), '{}'.format(lon)
obs.elevation = alt

# speed of light
c = 299792458.0 # m/s

# This sample rate is used because it is a multiple of the baud rate
# Also it allows capturing the whole 1 MHz channel, which ensures getting
# both orbcomm channels without tuning to a different center frequency
# However, it is a much higher sample rate than would be needed.
sample_rate = 1.2288e6
center_freq = 137.5e6

# dictionary/list to hold satellite plots
sat_gps_dict = {}
sat_plot_lines = []

# receive samples that are an integer multiple of 1024 from the RTLSDR
should_finish = False
queue_max_size = 30

# UDP ip address and port to listen on for samples from rtl_sdr
UDP_IP = ip_address
UDP_PORT = port_number

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM, socket.IPPROTO_UDP)
sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
sock.bind((UDP_IP, UDP_PORT))

# This function collects samples from the remote computer
def socket_get_samples(queue, context):

    sat = context['sat']
    sat_name = context['sat_name']
    obs = context['observer']

    while 1:
        try:
            # rtl_sdr over netcat sends small batches of samples at a time which
            # fills up the buffers to quickly. Concatenating many batches into
            # one bigger buffer allows faster processing.
            samples = np.array([], dtype=np.complex64)
            for xx in range(20):
                data, addr = sock.recvfrom(65565)
                # There are some very small packets that don't have samples in them
                if len(data) < 100:
                    continue

                # The samples are a string of binary data, uint8's, interleaved IQ
                # This code converts them to numpy complex64 samples
                temp_samples = np.fromstring(data, dtype=np.uint8)
                temp_samples = temp_samples.astype(np.float64)/127.5
                temp_samples = temp_samples.view(np.complex128) - (1 + 1j)
                samples = np.concatenate((samples, temp_samples))

            # Calculate the doppler shift for the satellite
            obs.date = ephem.now()
            sat.compute(obs)
            relative_vel = sat.range_velocity
            doppler = c/(c+relative_vel) * sat_center_frequency - sat_center_frequency

            # This code will catch the case when the user closes the matplotlib figure
            if queue.full():
                # empty out the queue so that the threads can join without error
                while queue.qsize() > 0:
                    queue.get()
                queue.put((None, None))
                break

            # If the satellite goes below the horizon, stops receiving samples
            if degrees(sat.alt) < min_elevation:
                queue.put((None, None))
                break

            # I pass a new contex to the processing thread
            new_context = {
                    'doppler':doppler,
                    'sat_name':sat_name,
                    'elevation':degrees(sat.alt),
                    'azimuth':degrees(sat.az)
                    }
            queue.put((samples, new_context))

        except KeyboardInterrupt:
            break
        # end while
    return 0

# This function runs in a separate process
# It collects samples from a queue, processes them and plots the results
def process_samples(queue):
    global decoder
    global should_finish
    global sat_gps_dict
    global sat_plot_lines

    ##########################################
    # configure the matplotlib plots
    if realtime_plotting:
        # Create a figure with 4 subplots
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
        img = plt.imread("./map.jpg")
        img_plot = ax_arr[1,1].imshow(img, extent=[-180, 180, -90, 90])
        loc_plot, = ax_arr[1,1].plot(lon, lat, 'ro')

        ax_arr[0,1].axis('off')
        plt.pause(0.0001)
    ##########################################

    tic = time()
    first_time = True

    while 1:
        try:
            samples, context = queue.get()
            if context == None:
                break

            doppler = context['doppler']
            sat_name = context['sat_name']
            az = context['azimuth']
            el = context['elevation']

            if first_time:
                decoder.first_samples(samples, doppler=doppler)
                first_time = False
                continue

            tic = time()
            packet_list = decoder.decode_samples_to_packets(samples, doppler)
            # if should_save_to_file:
            #     do saving stuff
            decoder.parse_packets()

            if realtime_plotting:
                if len(decoder.decimated_samples) == 0:
                    # This happens when the receiver resets
                    continue

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
                            sat_plot_lines.append(ax_arr[1,1].plot([0,], 'D:', color='cyan', markevery=[-1])[0])
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
                time_ratio = sample_rate/len(samples) * (toc-tic)

                bad_percent = decoder.bad_packets / (decoder.bad_packets + decoder.good_packets + 0.1)
                good_percent = decoder.good_packets / (decoder.bad_packets + decoder.good_packets + 0.1)

                text0 = ax_arr[0,1].text(0.10, 0.90, "Satellite Name: {}".format(sat_name))
                text1 = ax_arr[0,1].text(0.10, 0.80, "Lat/Lon: {:6.3f} {:6.3f}".format(decoder.sat_lon, decoder.sat_lat))
                text2 = ax_arr[0,1].text(0.10, 0.70, "Azimuth: {:6.0f}   Elevation: {:6.0f}".format(az, el))
                text3 = ax_arr[0,1].text(0.10, 0.60, "Bad Packets %:  {:4.1f}".format(bad_percent*100.0))
                text4 = ax_arr[0,1].text(0.10, 0.50, "Good Packets %: {:4.1f}".format(good_percent*100.0))
                text5 = ax_arr[0,1].text(0.10, 0.40, "Time ratio:     {:4.2f}".format(time_ratio))

                ax_arr[0,1].draw_artist(ax_arr[0,1].patch)
                ax_arr[0,1].draw_artist(text0)
                ax_arr[0,1].draw_artist(text1)
                ax_arr[0,1].draw_artist(text2)
                ax_arr[0,1].draw_artist(text3)
                ax_arr[0,1].draw_artist(text4)
                ax_arr[0,1].draw_artist(text5)

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
            print("Receiving from: {}".format(sat_name))
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

            print("Create multiprocessing queue")
            queue = mp.Queue(queue_max_size)

            context_dict = {
                            'observer':obs,
                            'sat':sat,
                            'sat_name':sat_name,
                            }

            print("Creating socket receiver process")
            s = mp.Process(target=socket_get_samples, args=(queue, context_dict))
            s.start()

            result = process_samples(queue)

            # After processing finishes, allow the socket process to terminate
            while s.is_alive():
                sleep(0.1)
            queue.close()
            queue.join_thread()
            s.join()

            # This is true if the user closed the plotting window.
            # It is not true if the satellite went below the horizon.
            if should_finish:
                break

        else:
            # If no satellite is overhead, find the next one that will be
            sat_detected = False

            # lon += 10.0
            # obs.lon = '{}'.format(lon)
            # continue

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

print("Program exited.")