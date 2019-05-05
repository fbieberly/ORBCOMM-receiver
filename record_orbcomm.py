
from time import time, sleep
from math import floor, degrees

import ephem
import numpy as np
from rtlsdr import RtlSdr
from scipy.io import savemat

from helpers import get_tle_lines
from sat_db import active_orbcomm_satellites


# Create a pyephem sat object for all the active satellites
# using latest TLE data
for name in active_orbcomm_satellites:
    sat_line0, sat_line1, sat_line2 = get_tle_lines(name, tle_dir='./tles')
    sat = ephem.readtle(sat_line0, sat_line1, sat_line2)
    active_orbcomm_satellites[name]['sat_obj'] = sat
    active_orbcomm_satellites[name]['tles'] = [sat_line0, sat_line1, sat_line2]


# PyEphem observer 
# Imput your receivers latitude, longitude and altitude
lat = 43.802953
lon = -99.210731
alt = 0 
obs = ephem.Observer()
obs.lat, obs.lon = '{}'.format(lat), '{}'.format(lon)
obs.elevation = alt # Technically is the altitude of observer
min_elevation = 30.0 # in degrees, in urban or heavily wooded areas, increase as appropriate

# Set RTLSDR parameters and initialize
sample_rate = 1.2288e6
center_freq = 137.5e6
gain = 'auto' # Use AGC

sdr = RtlSdr()
sdr.rs = sample_rate
sdr.gain = gain
sdr.fc = center_freq

# Set recording parameters
record_length = 2.0 # seconds
record_interval = 5.0 # seconds

# receive samples that are an integer multiple of 1024 from the RTLSDR
num_samples_per_recording = int((floor(record_length*sample_rate/1024)+1)*1024)

file_count = 100
while 1:
    try:
        start_loop = time()
        obs.date = ephem.now()

        sats = []
        tles = []
        for sat_name in active_orbcomm_satellites:
            sat = active_orbcomm_satellites[sat_name]['sat_obj']
            sat.compute(obs)
            if degrees(sat.alt) > min_elevation:
                sats.append(sat_name)
                tles.append(active_orbcomm_satellites[sat_name]['tles'])

        if len(sats) > 0:
            print("\nSatellites overhead: ")
            for sat_name in sats:
                sat = active_orbcomm_satellites[sat_name]['sat_obj']
                sat.compute(obs)
                print('{:20}: {:3.1f} degrees elevation'.format(sat_name, degrees(sat.alt)))
            record_time = time()
            samples = sdr.read_bytes(num_samples_per_recording*2)
            complex_samples = np.frombuffer(samples, dtype=np.uint8)
            complex_samples = complex_samples.astype(np.float32) - 128
            complex_samples = complex_samples.astype(np.float32).view(np.complex64)
            filename = '{:.3f}'.format(record_time).replace('.', 'p') + '.mat'
            save_dict = {
                        'samples':complex_samples,
                        'timestamp':record_time,
                        'sats': sats,
                        'tles': tles,
                        'fs': sample_rate,
                        'fc': center_freq,
                        'lat':lat,
                        'lon':lon,
                        'alt':alt,
                        }
            savemat('./data/' + filename, save_dict, do_compression=True)
            print("File saved: {}".format('./data/' + filename))
            file_count -= 1
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
                            print("Minutes until next satellite visible: {:.0f}".format(minute))
                            sleep(60)

                        break
                if sat_detected:
                    break
            if sat_detected == False:
                print("No upcoming satellite passes detected within 12 hours.")
                exit()


        if file_count <= 0:
            break
        loop_time_remaining = record_interval - (time()-start_loop)
        if loop_time_remaining > 0:
            sleep(loop_time_remaining)

    except KeyboardInterrupt:
        break

sdr.close()