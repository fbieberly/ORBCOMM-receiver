
from glob import glob

import numpy as np
from scipy.signal import butter, lfilter


def get_tle_lines(sat_name, tle_dir='./tles'):
    line0, line1, line2 = '', '', ''
    tles = glob(tle_dir + '/*.txt')
    for tle in tles:
        with open(tle, 'r') as f:
            text = 'abc'
            while text != '':
                text = f.readline()
                if sat_name.lower() in text.lower():
                    line0 = text
                    line1 = f.readline()
                    line2 = f.readline()
                    break
        if line0 != '':
            break
    return line0, line1, line2


def quad_interp(arr):
    k = np.argmax(arr)
    y1 =  abs(arr[k - 1])
    y2 =  abs(arr[k])
    y3 =  abs(arr[k + 1])
    d  = (y3 - y1) / (2 * (2 * y2 - y1 - y3))
    interp_k =  k + d
    return d

def butter_lowpass_filter(data, cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = lfilter(b, a, data)
    return y

def complex_mix(arr, freq_shift, sample_rate):
    duration = len(arr)*1.0/sample_rate
    t = np.arange(0, duration*2.0, 1.0/sample_rate)[:len(arr)] # Need to make sure you have as many of these samples as you have IQ samples. Sometimes I use: t = np.arange(0, duration*2, 1.0/sample_rate)[:len(signal)]
    complex_cos = np.exp( 1j * 2*np.pi * freq_shift * t, dtype=np.complex64)
    shifted_signal = arr * complex_cos
    return shifted_signal