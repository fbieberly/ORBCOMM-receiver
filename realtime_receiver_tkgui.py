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
# There is a "GUI" composed of a tkinter window and two matplotlib plots.
# The GUI has a few slightly useful bits of information for troubleshooting
# if you aren't getting a good reception. 
# The code uses multiprocessing to run the decoder in a separate thread from
# the RTLSDR code. I haven't gotten around to getting the code to finish
# cleanly. 
# 
#
##############################################################################


tk_imported = False
try:
    import tkinter as tk
    tk_imported = True
except ImportError:
    pass
try:
    import Tkinter as tk
    tk_imported = True
except ImportError:
    pass
if tk_imported == False:
    print("No Tkinter library found. Exiting.")
    exit()

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

import cPickle as pickle
from time import time, sleep
import multiprocessing as mp
from math import floor, degrees, log10

import ephem
import numpy as np
from rtlsdr import RtlSdr
from scipy.io import savemat
from scipy.signal import welch
import matplotlib.pyplot as plt


from helpers import get_tle_lines
from realtime_decoder import RealtimeDecoder
from sat_db import active_orbcomm_satellites

#####
# Tkinter gui

window = tk.Tk()
window.title("ORBCOMM Decoder")
window.geometry('450x200')
sat_name_lbl = tk.Label(window, text="Satellite Name: ")
sat_name_lbl.grid(column=0, row=0, sticky=tk.W)

angle_lbl = tk.Label(window, text="Phase angle:  /")
angle_lbl.grid(column=0, row=1, sticky=tk.W)

good_pack_lbl = tk.Label(window, text="Good packets: ")
good_pack_lbl.grid(column=0, row=2, sticky=tk.W)

bad_pack_lbl = tk.Label(window, text="Bad packets: ")
bad_pack_lbl.grid(column=0, row=3, sticky=tk.W)

time_ratio_lbl = tk.Label(window, text="Time ratio: ")
time_ratio_lbl.grid(column=0, row=4, sticky=tk.W)
window.update()

######
# configure the matplotlib plots
const_fig, const_ax = plt.subplots()
cost_line, = const_ax.plot([0,], 'bx')
plt.xlim((-2, 2))
plt.ylim((-2, 2))

plt.pause(0.0001)
fft_fig, ftt_ax = plt.subplots()
fft_line, = ftt_ax.plot([0,], 'b')
plt.xlim((-6e3, 6e3))
plt.ylim((-50, 20))
plt.pause(0.0001)

# Create a pyephem sat object for all the active satellites
# using latest TLE data
for name in active_orbcomm_satellites:
    sat_line0, sat_line1, sat_line2 = get_tle_lines(name, tle_dir='./tles')
    sat = ephem.readtle(sat_line0, sat_line1, sat_line2)
    active_orbcomm_satellites[name]['sat_obj'] = sat

# PyEphem observer 
# Imput your receivers latitude, longitude and altitude
lat = 47.245255
lon = -122.469978
alt = 0 
obs = ephem.Observer()
obs.lat, obs.lon = '{}'.format(lat), '{}'.format(lon)
obs.elevation = alt
min_elevation = 5.0 # in degrees, in urban or heavily wooded areas, increase as appropriate

# speed of light
c = 299792458.0 # m/s

# Set RTLSDR parameters and initialize
sample_rate = 1.2288e6
center_freq = 137.5e6
gain = 'auto' # Use AGC

sdr = RtlSdr()
sdr.rs = sample_rate
sdr.gain = gain
sdr.fc = center_freq

# receive samples that are an integer multiple of 1024 from the RTLSDR
num_samples_per_recording = int(1024*256)

def rtlsdr_callback(samples, context):
    global decoder
    sdr = context['sdr']
    sat = context['sat']
    sat_name = context['sat_name']
    obs = context['observer']

    obs.date = ephem.now()
    if degrees(sat.alt) < min_elevation:
        queue.put((None, None))
        sdr.cancel_read_async()
        return 0

    obs.date = ephem.now()
    sat.compute(obs)
    relative_vel = sat.range_velocity
    doppler = c/(c+relative_vel) * sat_center_frequency - sat_center_frequency

    new_context = {
            'doppler':doppler,
            'sat_name':sat_name,
            }

    queue.put((samples, new_context))

def process_samples(queue):
    global decoder

    while 1:
        samples, context = queue.get()
        if context == None:
            break

        doppler = context['doppler']
        sat_name = context['sat_name']

        tic = time()
        decoder.decode_samples_to_packets(samples, doppler)
        decoder.parse_packets()

        f, pxx = welch(decoder.decimated_samples, fs=decoder.sample_rate/decoder.decimation, \
                          return_onesided=False, nperseg= int(len(decoder.decimated_samples)/10), scaling='density')
        f = (np.roll(f, int(len(f)/2)))
        pxx = np.roll(pxx, int(len(pxx)/2))
        pxx = 20*np.log10(pxx)
        pxx -= np.median(pxx)

        symbols = decoder.symbols
        symbols /= np.median(np.abs(symbols))

        cost_line.set_xdata(symbols.real)
        cost_line.set_ydata(symbols.imag)
        const_ax.draw_artist(const_ax.patch)
        const_ax.draw_artist(cost_line)
        const_fig.canvas.update()
        const_fig.canvas.flush_events()

        fft_line.set_xdata(f)
        fft_line.set_ydata(pxx)
        ftt_ax.draw_artist(ftt_ax.patch)
        ftt_ax.draw_artist(fft_line)
        fft_fig.canvas.update()
        fft_fig.canvas.flush_events()

        toc = time()
        time_ratio = sdr.sample_rate/len(samples) * (toc-tic)

        try:
            bad_percent = decoder.bad_packets / (decoder.bad_packets + decoder.good_packets + 0.1)
            good_percent = decoder.good_packets / (decoder.bad_packets + decoder.good_packets + 0.1)
            ave_above = decoder.ave_angles_above_zero
            ave_below = decoder.ave_angles_below_zero
            sat_name_lbl.configure(   text="Satellite Name: {}".format(sat_name), justify=tk.LEFT)
            angle_lbl.configure(      text="Phase angle:    {:} / {:}".format(int(ave_above), int(ave_below)), justify=tk.LEFT)
            good_pack_lbl.configure(  text="Good packets:   {:}%".format(int(good_percent*100)), justify=tk.LEFT)
            bad_pack_lbl.configure(   text="Bad packets:    {:}%".format(int(bad_percent*100)), justify=tk.LEFT)
            time_ratio_lbl.configure( text="Time ratio:     {:4.2f}".format(time_ratio), justify=tk.LEFT)

            window.update_idletasks()
            window.update()
        except tk.TclError:
            break

while 1:
    try:
        start_loop = time()
        obs.date = ephem.now()

        sats = []
        for sat_name in active_orbcomm_satellites:
            sat = active_orbcomm_satellites[sat_name]['sat_obj']
            sat.compute(obs)
            if degrees(sat.alt) > min_elevation:
                sats.append((sat_name, sat, degrees(sat.alt)))

        if len(sats) > 0:
            print("\nSatellites overhead: ")
            sorted_sats = sorted(sats, key=lambda x: x[2], reverse=True)
            for sat_name, sat, degrees_above_horizon in sorted_sats:
                print('{:20}: {:3.1f} degrees elevation'.format(sat_name, degrees_above_horizon))

            sat_name = sorted_sats[0][0]
            sat = sorted_sats[0][1]
            frequencies = active_orbcomm_satellites[sorted_sats[0][0]]['frequencies']
            print("Satellite frequencies: {}".format(frequencies))
            sat_center_frequency = frequencies[1]
            
            # Decoding the lower of the two frequencies
            decoder = RealtimeDecoder(sat_center_frequency, center_freq=center_freq, sample_rate=sample_rate)


            obs.date = ephem.now()
            sat.compute(obs)
            relative_vel = sat.range_velocity
            doppler = c/(c+relative_vel) * sat_center_frequency - sat_center_frequency

            print('Recording samples.')
            # Record samples twice
            samples = sdr.read_samples(num_samples_per_recording)
            samples = sdr.read_samples(num_samples_per_recording)

            print("Create multiprocessing queue")
            # create multiprocessing queue
            queue = mp.Queue()

            print("Creating multiprocessing process")
            p = mp.Process(target=process_samples, args=(queue,))
            p.start()

            decoder.first_samples(samples, doppler=doppler)

            context_dict = {
                            'sdr':sdr,
                            'observer':obs,
                            'sat':sat,
                            'sat_name':sat_name,
                            }
            sdr.read_samples_async(rtlsdr_callback, num_samples_per_recording, context_dict)
            
            queue.close()
            queue.join_thread()
            p.join()

        else:
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