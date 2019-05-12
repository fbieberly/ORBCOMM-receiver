
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

# From: https://github.com/veeresht/CommPy/blob/master/commpy/filters.py
def rrcosfilter(N, alpha, Ts, Fs):
    """
    Generates a root raised cosine (RRC) filter (FIR) impulse response.
    
    Parameters
    ----------
    N : int 
        Length of the filter in samples.
    
    alpha: float
        Roll off factor (Valid values are [0, 1]).
    
    Ts : float
        Symbol period in seconds.
    
    Fs : float 
        Sampling Rate in Hz.
    
    Returns
    ---------
    h_rrc : 1-D ndarray of floats
        Impulse response of the root raised cosine filter.
    
    time_idx : 1-D ndarray of floats 
        Array containing the time indices, in seconds, for 
        the impulse response.
    """
    N = int(N)
    T_delta = 1/float(Fs)
    time_idx = ((np.arange(N)-N/2))*T_delta
    sample_num = np.arange(N)
    h_rrc = np.zeros(N, dtype=float)
        
    for x in sample_num:
        t = (x-N/2)*T_delta
        if t == 0.0:
            h_rrc[x] = 1.0 - alpha + (4*alpha/np.pi)
        elif alpha != 0 and t == Ts/(4*alpha):
            h_rrc[x] = (alpha/np.sqrt(2))*(((1+2/np.pi)* \
                    (np.sin(np.pi/(4*alpha)))) + ((1-2/np.pi)*(np.cos(np.pi/(4*alpha)))))
        elif alpha != 0 and t == -Ts/(4*alpha):
            h_rrc[x] = (alpha/np.sqrt(2))*(((1+2/np.pi)* \
                    (np.sin(np.pi/(4*alpha)))) + ((1-2/np.pi)*(np.cos(np.pi/(4*alpha)))))
        else:
            h_rrc[x] = (np.sin(np.pi*t*(1-alpha)/Ts) +  \
                    4*alpha*(t/Ts)*np.cos(np.pi*t*(1+alpha)/Ts))/ \
                    (np.pi*t*(1-(4*alpha*t/Ts)*(4*alpha*t/Ts))/Ts)
        
    return time_idx, h_rrc

def fletcher_checksum(hex_data_str):
    sum1 = 0
    sum2 = 0

    if len(hex_data_str)%2 == 1:
        # hex_data_str = '0' + hex_data_str
        hex_data_str += '0'
        print('...')

    for xx in range(0, len(hex_data_str)-1, 2):
        val = int(hex_data_str[xx:xx+2], 16)
        sum1 = (sum1 + val)%256
        sum2 = (sum1 + sum2)%256

    return '{:02X}{:02X}'.format(sum2, sum1)

def reverse_endian(hex_data_str):
    out_string = ''
    for xx in range(0, len(hex_data_str)-1, 2):
        val = '{:08b}'.format(int(hex_data_str[xx:xx+2], 16))
        out_string += '{:02X}'.format(int(val[::-1],2))

    return out_string

def ecef_to_lla(x_ecef, y_ecef, z_ecef):
    # From: http://www.epsg.org/Portals/0/373-07-2.pdf?ver=2019-03-08-165437-017
    # page 97
    f = 1.0 / 298.257223563

    wgs84_a = 6378137.0 # m
    wgs84_b = wgs84_a*(1.-f)
    wgs84_e_sqrd = 1. - (wgs84_b**2)/(wgs84_a**2)

    wgs84_eps = wgs84_e_sqrd / (1 - wgs84_e_sqrd)
    wgs84_p = np.sqrt(x_ecef**2 + y_ecef**2)
    wgs84_q = np.arctan2((z_ecef * wgs84_a), (wgs84_p * wgs84_b))

    phi = np.arctan2((z_ecef + wgs84_eps * wgs84_b * np.sin(wgs84_q)**3), \
                    (wgs84_p - wgs84_e_sqrd * wgs84_a * np.cos(wgs84_q)**3))
    lamd = np.arctan2(y_ecef, x_ecef)

    wgs84_v = (wgs84_a / np.sqrt(1 - wgs84_e_sqrd * np.sin(phi)**2))
    h = (wgs84_p / np.cos(phi)) - wgs84_v

    return np.degrees(phi), np.degrees(lamd), h
