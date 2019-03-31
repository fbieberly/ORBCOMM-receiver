import numpy as np
import scipy.signal as scisig
import matplotlib.pyplot as plt

def gen_bits(num_bits, zero_runs=False, one_runs=False):
    bits = np.random.randint(0, 2, size=num_bits)
    # Can add runs of 10 consecutive zeros
    if zero_runs:
        for xx in range(1, int(len(bits)/30)):
            bits[xx * 30 - 20:xx * 30 - 10] = np.zeros(10)
    # Can add runs of 10 consecutive ones
    if one_runs:
        for xx in range(1, len(bits)/30):
            bits[xx * 30 - 10:xx * 30] = np.ones(10)
    return bits

def bpsk_symbols(bits, samples_per_symbol):
    """
    inputs:
    bits: a numpy array of 0, 1 integers
    samples_per_symbol: an integer
    output:
    symbols: a numpy array of np.complex64 samples
    """
    symbols = np.zeros(int(len(bits) * samples_per_symbol), dtype=np.complex64)
    symbols.real[::samples_per_symbol] = (bits * 2.0) - 1
    return symbols

def qpsk_symbols(bits, samples_per_symbol):
    """
    inputs:
    bits: a numpy array of 0, 1 integers
    samples_per_symbol: an integer
    output:
    symbols: a numpy array of np.complex64 samples
             array is zero stuffed up to samples_per_symbol rate
    """
    symbols = np.zeros(int(len(bits)/2 * samples_per_symbol), dtype=np.complex64)
    symbols.real[::samples_per_symbol] = (bits[0::2] * 2.0) - 1
    symbols.imag[::samples_per_symbol] = (bits[1::2] * 2.0) - 1
    return symbols

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




if __name__ == '__main__':
    num_symbols = 500
    samples_per_symbol = 8

    bits = gen_bits(num_symbols)
    symbols = qpsk_symbols(bits, samples_per_symbol)

    alpha = 0.4
    baudrate = 1.
    num_of_symbols_half_filter = 4
    rrc_num_taps = samples_per_symbol * num_of_symbols_half_filter * 2.
    t_idx, rrc_taps = rrcosfilter(rrc_num_taps, alpha, baudrate, samples_per_symbol)

    plt.figure()
    plt.title("RRC taps")
    plt.plot(rrc_taps)
    plt.grid()

    # Normalize RRC filter
    rrc_taps /= np.sum(rrc_taps)

    # Pulse shape filter the tx_samples
    # Notice the rrc filter taps are multiplied by samples_per_symbol
    # This normalizes the power in the samples (after the zero stuffing)
    tx_samples = scisig.lfilter(rrc_taps[::]*samples_per_symbol, [1.0], symbols[::])

    # interpolate tx samples to have timing error
    td = 1.5
    tau = td                               # timing correction to be applied
    th  = np.array([-2, -1, 0, 1, 2])       # time vector for interpolating filter
    w   = np.sqrt(np.hamming(5).T)          # window function for interpolating filter
    hi  = np.sinc(th + tau) * w             # interpolating filter coefficients
    hi /= np.sum(hi)
    
    tx_samples = scisig.lfilter(hi, [1.], tx_samples)     # interpolate matched filter output

    # MF the samples at the receiver
    mf_samples = scisig.lfilter(rrc_taps[::], [1.0], tx_samples[::])

    # interpolate rx samples to remove timing error
    tau = -td                               # timing correction to be applied
    th  = np.array([-2, -1, 0, 1, 2])       # time vector for interpolating filter
    w   = np.sqrt(np.hamming(5).T)          # window function for interpolating filter
    hi  = np.sinc(th + tau) * w             # interpolating filter coefficients
    hi /= np.sum(hi)
    mf_samples_timing = scisig.lfilter(hi, [1.], mf_samples)     # interpolate matched filter output


    plt.figure()
    plt.subplot(121)
    plt.title("QPSK symbols, RRC filtering")
    plt.xlabel("Samples")
    plt.ylabel("Magnitude")
    plt.scatter(tx_samples.real[2::samples_per_symbol], \
                tx_samples.imag[2::samples_per_symbol], marker='x', label='Pulse shaped')

    plt.scatter(mf_samples.real[2::samples_per_symbol], \
                mf_samples.imag[2::samples_per_symbol], marker='x', label='Match filtered')

    plt.scatter(mf_samples_timing.real[4::samples_per_symbol], \
                mf_samples_timing.imag[4::samples_per_symbol], marker='x', label='Timing recovery')

    plt.scatter(symbols.real[::samples_per_symbol], \
                symbols.imag[::samples_per_symbol], marker='x', label='Symbols')
    plt.legend(loc='best')


    plt.subplot(122)
    plt.title("QPSK symbols, After MF")
    plt.xlabel("I")
    plt.ylabel("Q")
    plt.plot(mf_samples.real[:250])
    plt.plot(mf_samples.imag[:250])
    plt.show()