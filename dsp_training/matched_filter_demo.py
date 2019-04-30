##############################################################################
#
# Author: Frank Bieberly
# Date: 22 March 2019
# Name: matched_filter_demo.py
# Description: 
# This is a demo of transmitting and receiving qpsk samples. It uses a RRC 
# pulse shaping filter and matched filter on the receiver. Note that the 
# Tx signal has ISI after the RRC filter, but matched filtering at the receiver
# removes almost all of the ISI.
#
##############################################################################

import numpy as np
import scipy.signal as scisig
import matplotlib.pyplot as plt
from utilities import gen_bits, qpsk_symbols, rrcosfilter


if __name__ == '__main__':
    num_symbols = 500
    samples_per_symbol = 4

    bits = gen_bits(num_symbols, zero_runs=False, one_runs=False)
    symbols = qpsk_symbols(bits, samples_per_symbol)

    # Create RRC taps
    alpha = 0.4
    baudrate = 1.
    num_of_symbols_half_filter = 4
    rrc_num_taps = samples_per_symbol * num_of_symbols_half_filter * 2.
    t_idx, rrc_taps = rrcosfilter(rrc_num_taps, alpha, baudrate, samples_per_symbol)

    # Plot the RRC taps
    plt.figure()
    plt.title("RRC taps")
    plt.plot(rrc_taps)
    plt.grid()

    # Normalize RRC filter, helps with the plotting
    rrc_taps /= np.sum(rrc_taps)

    # Pulse shape filter the tx_samples
    # Notice the rrc filter taps are multiplied by samples_per_symbol
    # This normalizes the power in the samples (after the zero stuffing)
    tx_samples = scisig.lfilter(rrc_taps[::]*samples_per_symbol, [1.0], symbols[::])

    # Matched filter the samples at the receiver with same RRC filter
    mf_samples = scisig.lfilter(rrc_taps[::], [1.0], tx_samples[::])

    # Plot the symbols, Tx symbols, and Rx symbols
    plt.figure()
    plt.subplot(121)
    plt.title("QPSK symbols, RRC filtering")
    plt.xlabel("Samples")
    plt.ylabel("Magnitude")
    plt.scatter(tx_samples.real[::samples_per_symbol], \
                tx_samples.imag[::samples_per_symbol], marker='x', label='Tx after filter')

    plt.scatter(mf_samples.real[::samples_per_symbol], \
                mf_samples.imag[::samples_per_symbol], marker='x', label='Rx after filter')

    plt.scatter(symbols.real[::samples_per_symbol], \
                symbols.imag[::samples_per_symbol], marker='x', label='Symbols')
    plt.legend(loc='best')

    # Plot the time domain signal
    plt.subplot(122)
    plt.title("QPSK symbols, After MF")
    plt.xlabel("I")
    plt.ylabel("Q")
    plt.plot(mf_samples.real[:250])
    plt.plot(mf_samples.imag[:250])
    plt.show()