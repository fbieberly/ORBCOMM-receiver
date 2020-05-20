# ORBCOMM receiver
A software receiver for ORBCOMM satellite transmissions.  

Please read the [wiki] for more information.

[wiki]: https://github.com/fbieberly/ORBCOMM-receiver/wiki  

## Description

This is a software receiver for decoding packets from ORBCOMM satellites. For known packet types the packet data is decoded.  

I am writing this decoder as an instructional personal project. Hopefully it
can be used by others to learn about designing satellite communication
receivers.  

If you want a more full-featured ORBCOMM receiver please check out:  
https://www.coaa.co.uk/orbcommplotter.htm  
http://f6cte.free.fr/index_anglais.htm  



## Dependencies

Should work with either Python 2.X or 3.X

I use:  
[pyrtlsdr] to record the RF signal with an RTLSDR receiver.  
[NumPy] and [SciPy] are used for signal processing.  
[PyEphem] is used to calculate Az/El and doppler shift of the satellites.  
[Matplotlib] is used for plotting data.

pip install pyrtlsdr numpy scipy pyephem matplotlib



[PyEphem]: https://rhodesmill.org/pyephem/index.html
[NumPy]: https://www.numpy.org/
[SciPy]: https://www.scipy.org/
[pyrtlsdr]: https://github.com/roger-/pyrtlsdr
[Matplotlib]: https://matplotlib.org/

## Getting started

#### Offline recording and decoding
1. First run the _update_orbcomm_tle.py_ script to get the latest two-line elements for the ORBCOMM satellites.
2. Update latitude and longitude of your receiver in _CONFIG.py_
    1. Should be in the format: ```lat, lon = 40.919295, -74.087683```
3. Record IQ data by running _record_orbcomm.py_
    1. Alternatively, use the _record_orbcomm_long.py_ script which will produce a longer recording that has been decimated to save storage space.
    1. The longest recording that _record_orbcomm.py_ can make is about 3 seconds, _record_orbcomm_long.py_ can make recordings that are ~100 seconds long.
4. Run _file_decoder.py_ to decode a single recording file (defaults to the first file in the /data folder)
    1. The file it decodes is selected near the top of the file (Line 37). Change it there if you wish to decode other files.
    1. Note: There are 2 data files that come with this repo. If you want to process your own recordings, you'll need to remove those files (or just point the file_decoder.py script to your files).

```
EXAMPLE OUTPUT
Filename: ./data/1552071892p6.mat
Timestamp: 1552071892.6
Data collected on: 2019-03-08 19:04:52.600117
Satellites in recording: orbcomm fm114
SDR Sample rate: 1228800.0 Hz
SDR Center frequency: 137500000.0 Hz
Satellite frequencies: 137287500.0, 137737500.0
Number of samples: 2458624
Remaining frequency offset after doppler compensation: -141.0 Hz
Number of possible packets: 100
Bit stream offset: 3

List of packets: (### indicates checksum failed)
### Unrecognized packet: AF0C3958A56A1A7A9BE0ED2B
Fill: data: 19F1D8528E1EF9701DED 
Fill: data: 5A8C1A5E354CE775C6A3 
Fill: data: 6FE4A5F1DDC8B12CADE5 
Fill: data: 567DC3683187453D728D 
Fill: data: 463E5D5FB36A15E68001 
Message: msg_packet_num: 0 msg_total_length: 2 data: 001F05CE01C0721828 
Message: msg_packet_num: 1 msg_total_length: 2 data: 507102000000000032 
Message: msg_packet_num: 0 msg_total_length: 2 data: 1167036942905C869C 
Message: msg_packet_num: 1 msg_total_length: 2 data: C10183D10BE0A32C54 
Message: msg_packet_num: 0 msg_total_length: 3 data: A241000129687B035E 
Message: msg_packet_num: 1 msg_total_length: 3 data: 921E026637E0228277 
Message: msg_packet_num: 2 msg_total_length: 3 data: A236830000000000F7 
Message: msg_packet_num: 0 msg_total_length: 2 data: 836B01E54270B20226 
Message: msg_packet_num: 1 msg_total_length: 2 data: 734502B42B00000048 
Fill: data: B3489015957F14F39CC3 
Fill: data: D9B2B6A2EB952C9A23AD 
Unrecognized packet: 0B01FD24CCCCCC204501CF3A
Sync: code: 65A8F9 sat_id: 2C 
Downlink_info: msg_packet_num: 0 msg_total_length: 3 data: 27310750A005640094 
Downlink_info: msg_packet_num: 1 msg_total_length: 3 data: 7D000BB89010130195 
Downlink_info: msg_packet_num: 2 msg_total_length: 3 data: 1D011400000000003F 
Network: msg_packet_num: 0 msg_total_length: 1 data: 7800010000000000E2 
Ephemeris: sat_id: 2C data: 98E3D5043B9BC34CDDF04C3CE66F98D5A307FB07E1D8 
    Current satellite time: 2019-03-08 19:04:53 Z
    Lat/Lon:        44.9555, -116.0878, Altitude:  715.3 km, Velocity: 7151.4 m/s
    Ephem Lat/Lon:  44.7939, -116.0389, Altitude:  715.2 km, Velocity: 7170.9 m/s
    Difference in reported and ephemeris position:   21.0 km
    Difference in reported and ephemeris velocity:   19.5 m/s
Unrecognized packet: 0B011C2AEDEEEE409B02A761
Unrecognized packet: 0B01641D76989944450169D9
Unrecognized packet: 0A48B99524C221A30012BCE8
Unrecognized packet: 0A4C3F01207164CF01127E15
Unrecognized packet: 0A506D28227274970012FF61
Unrecognized packet: 0A54E4BB7A6A5499021662B8
```



#### Real-time recording and decoding
1. First run the _update_orbcomm_tle.py_ script to get the latest two-line elements for the ORBCOMM satellites.
2. Update latitude and longitude of your receiver in _CONFIG.py_
    1. Should be in the format: ```lat, lon = 40.919295, -74.087683```
3. Run _realtime_receiver.py_
    1. If there is no satellite overhead it will tell you how long the wait is.
    1. I recommend you use [gPredict] to know where the ORBCOMM satellites are.
    1. Note that not all the ORBCOMM satellites still transmit. Look in _sat_db.py_ to see the active ones.  
  
  
#### Real-time recording and decoding from remote RTL-SDR
1. First run the _update_orbcomm_tle.py_ script to get the latest two-line elements for the ORBCOMM satellites.
1. Update latitude and longitude of your receiver in _CONFIG.py_
    1. Should be in the format: ```lat, lon = 40.919295, -74.087683```
1. Update ip_address and port_number of your 'processing computer' in CONFIG.py
    1. Should be on the same network as your remote computer.
1. Run _realtime_receiver_network.py_
    1. If there is no satellite overhead it will tell you how long the wait is.
    1. The program will block waiting for data to come in to the port.
1. Execute the following command on your remote computer (with the RTL-SDR): 
    1. ```rtl_sdr -f 137500000 -s 1228800 -g 0 - | netcat -uv ip_address port_number```
    1. ip_address and port_number should be the same as what you put in the CONFIG.py file.
    1. If latency is high to the remote computer, the program may not work well.

[gPredict]: http://gpredict.oz9aec.net/


## DSP Training

In the dsp_training folder are a number of scripts that I used to help me understand the DSP that I needed to decode the ORBCOMM signals. The scripts are simulation only and help understand phase recovery, timing recovery, creating symbols from bits, mixing, filtering, etc.



## Scripts


Scripts include:
- _sat_db.py_: just a dictionary of ORBCOMM satellites I know are active
- _orbcomm_packet.py_: a list of all the known ORBCOMM packet types and their components
- _helpers.py_: a file with useful helper functions
- _plot_recording_waterfall.py_: plots a waterfall of recordings
- _update_orbcomm_tle.py_: downloads the latest ORBCOMM tles from celestrack.com
- _record_orbcomm.py_: records ORBCOMM satellites when they are overhead with an RTLSDR
- record_orbcomm_long.py: same as record_orbcomm.py, but it decimates the signal so that it can make much longer recordings (only records one orbcomm channel at a time).
- _file_decoder.py_: If you have .mat files in the data folder, this script will attempt to decode one
- _realtime_decoder.py_: This is a class for doing decoding of a realtime stream of samples
- _realtime_receiver.py_: This is a script that does realtime decoding of the ORBCOMM signal, plus some interesting plots
- _realtime_receiver_network.py_: This is a script that does realtime decoding of the ORBCOMM signal, it receives samples from an RTL-SDR dongle that is on a remote computer.





## References

I used these two resources as my primary references.

http://mdkenny.customer.netspace.net.au/Orbcomm.pdf  
http://www.decodesystems.com/orbcomm.html  


## Data format

In the data folder is a couple files of samples that I have recorded.

The files are .mat files. They can be opened with MATLAB or Python (using SciPy's loadmat function).

The files include metadata:
- fc: center frequency
- fs: sample rate
- sats: a list of the names of the satellites overhead
- tles: a list of lists of the tle lines for each satellite (in the order of the sats list)
- timestamp: unix time of the start of the recording
- samples: a numpy complex64 array of the samples
- lat: the latitude of the receiver when the samples were recorded
- lon: the longitude of the receiver when the samples were recorded
- alt: the elevation of the receiver when the samples were recorded

Look at the file_decoder.py script to see an example of how to access the metadata.
