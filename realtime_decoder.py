##############################################################################
#
# Author: Frank Bieberly
# Date: 14 July 2019
# Name: realtime_decoder.py
# Description:
# This is a class designed for real time decoding of orbcomm satellite passes
#
#
##############################################################################

import numpy as np
from scipy.signal import firwin, welch
from datetime import datetime, timedelta

from orbcomm_packet import packet_dict
from helpers import complex_mix, rrcosfilter
from helpers import fletcher_checksum, ecef_to_lla


class RealtimeDecoder():

    def __init__(self, sat_freq, center_freq=137.5e6, sample_rate=1.2288e6):
        self.sample_rate = sample_rate
        self.center_frequency = center_freq
        self.sat_center_frequency = sat_freq
        self.init_default_values()

    def init_default_values(self):
        # Initialize all required variables
        # Some are hardcoded
        self.doppler = 0.0
        self.frequency_offset = 0.0
        self.symbol_offset = -1
        self.first_median = 1.0
        self.mix_phase = 0.0

        self.size_of_packets = int(12*8)
        self.bit_offset = -1
        self.baud_rate = 4800.0
        self.samples_per_symbol = 2
        self.decimation = int(self.sample_rate/(self.samples_per_symbol*self.baud_rate))

        self.bit_string = ''
        self.packets =[]
        self.reverse = False
        self.err = 0.0

        # Debug info
        self.ave_angles_above_zero = 0.0
        self.ave_angles_below_zero = 0.0
        self.good_packets = 0.0
        self.bad_packets = 0.0

        self.packet_dict = packet_dict
        self.packet_headers = [packet_dict[packet_type]['header'] for packet_type in packet_dict]

        # Timing recovery loop initialization
        self.tr_tau     = 0.                        # initial timing offset estimate
        self.tr_dtau    = 0.                        # initial timing _rate_ offset estimate
        self.tr_buf     = np.zeros(5, dtype=np.complex64) # filter buffer
        self.tr_th      = np.array([-2., -1., 0., 1., 2.]) # time vector for interpolating filter
        self.tr_w       = np.sqrt(np.hamming(5).T)  # window function for interpolating filter
        self.tr_alpha   = 0.0005                    # loop filter bandwidth (timing _rate_ adjustment factor)
        self.tr_beta    = 2*np.sqrt(self.tr_alpha)          # (timing _phase_ adjustment factor)
        self.tr_counter = 0         # interpolating filter adds delay
        self.tr_buf_dz  = [0.,0.,0.]
        self.tr_last_sample = (0 + 0j)
        self.last_symbol = (0 + 0j)

        # Initialize the costas loop parameters
        self.cl_alpha = 0.05
        self.cl_beta = 0.2 * self.cl_alpha**2
        self.cl_frequency_out = 0.0
        self.cl_phase_est = 0.0

        # Initialize LPF parameters
        self.lpf_filter_freq = 10e3
        self.lpf_order = 181   # filter order
        self.lpf_nyq = 0.5 * self.sample_rate
        self.lpf_normal_cutoff = self.lpf_filter_freq / self.lpf_nyq
        self.lpf_fir_taps = firwin(self.lpf_order, self.lpf_filter_freq, fs=self.sample_rate)

        # Create RRC taps
        rrc_alpha = 0.4
        rrc_baudrate = 1.
        rrc_num_of_symbols_half_filter = 8.
        self.rrc_num_taps = self.samples_per_symbol * rrc_num_of_symbols_half_filter * 2. + 1.
        rrc_t_idx, self.rrc_taps = rrcosfilter(self.rrc_num_taps, rrc_alpha, rrc_baudrate, self.samples_per_symbol)

        # sample buffers
        self.lpf_remaining_samples = np.array([])
        self.rrc_remaining_samples = np.array([])

        # store satellite info from packets
        self.sat_id = ''
        self.sat_lat = 0.0
        self.sat_lon = 0.0
        self.sat_alt = 0.0

        self.decimated_samples = np.array([])
        self.symbols = np.array([])
        self.mixed_down_samples = np.array([])


    def clean(self):
        # reset the decoder to receive a new satellite pass
        self.init_default_values()

    def first_samples(self, samples, doppler=0.0):
        # Use the first block of samples to estimate frequency offsets
        self.doppler=doppler
        self.first_median = np.median(np.abs(samples))
        norm_sample_block = samples / self.first_median

        # Mix the samples to baseband (compensating for doppler shift)
        # There will be a residual carrier error because of the RLTSDR frequency offset
        freq_shift = self.center_frequency - self.sat_center_frequency - self.doppler
        mixed_down_samples, _ = complex_mix(norm_sample_block, freq_shift, self.sample_rate)

        samps_to_filter = int(len(mixed_down_samples) - self.lpf_order)
        while samps_to_filter % self.decimation != 0:
            samps_to_filter -= 1

        # Low pass filter and decimate in one step
        decimated_samples = np.zeros(int(samps_to_filter/self.decimation), dtype=np.complex64)
        for yy in range(0, len(decimated_samples)):
            temp_samples = mixed_down_samples[yy*self.decimation:int(yy*self.decimation+self.lpf_order)]
            decimated_samples[yy] = np.dot(self.lpf_fir_taps, temp_samples)

        # estimate remaining carrier error (RTLSDR frequency error)
        # signal to the fourth power, take the fft, peak is at frequency offset
        rbw = 1
        nperseg = int(len(decimated_samples))
        signal_to_4th_power = np.power(decimated_samples, 4)
        f, pxx = welch(signal_to_4th_power, fs=self.sample_rate/self.decimation, \
                              return_onesided=False, nperseg=nperseg, scaling='density')
        f = (np.roll(f, int(len(f)/2)))
        pxx = np.roll(pxx, int(len(pxx)/2))
        search_window = int(1000. / ((self.sample_rate/self.decimation)/nperseg))  # search +/- 1 kHz around fc
        frequency_peak = np.argmax(pxx[int(len(pxx)/2 - search_window):int(len(pxx)/2 + search_window)])
        self.frequency_offset = -(frequency_peak - search_window)*(self.sample_rate/self.decimation/nperseg) / 4
        print("Frequency offset: {}".format(self.frequency_offset))

    def decode_samples_to_packets(self, samples, doppler=0.0):

        # on the very first call, first_samples will try to determine the frequency offset
        # of the SDR.
        if self.frequency_offset == None:
            self.first_samples(samples, doppler=doppler)

        # process samples
        # if doppler shift is provided, it is compensated for
        # otherwise, the phase recovery algorithm attempts to control it
        self.doppler = doppler
        freq_shift = self.center_frequency - self.sat_center_frequency - self.doppler + self.frequency_offset
        # norm_sample_block = samples / self.first_median
        norm_sample_block = samples / np.median(np.abs(samples))

        mixed_down_samples, self.mix_phase = complex_mix(norm_sample_block, freq_shift, self.sample_rate, self.mix_phase)
        mix_sample_buffer = np.concatenate([self.lpf_remaining_samples, mixed_down_samples])
        self.mixed_down_samples = mixed_down_samples

        samps_to_filter = int(len(mix_sample_buffer) - self.lpf_order)
        while samps_to_filter % self.decimation != 0:
            samps_to_filter -= 1

        # Low pass filter and decimate in one step
        self.decimated_samples = np.zeros(int(samps_to_filter / self.decimation), dtype=np.complex64)
        for yy in range(0, len(self.decimated_samples)):
            temp_samples = mix_sample_buffer[yy*self.decimation:int(yy * self.decimation + self.lpf_order)]
            self.decimated_samples[yy] = np.dot(self.lpf_fir_taps, temp_samples)
        self.lpf_remaining_samples = mix_sample_buffer[int((yy + 1) * self.decimation):]

        decim_sample_buffer = np.concatenate([self.rrc_remaining_samples, self.decimated_samples])
        matched_filtered_samples = np.zeros(int(len(decim_sample_buffer) - self.rrc_num_taps), dtype=np.complex64)
        for yy in range(len(matched_filtered_samples)):
            matched_filtered_samples[yy] = np.dot(self.rrc_taps, decim_sample_buffer[yy:int(yy+self.rrc_num_taps)])

        self.rrc_remaining_samples = decim_sample_buffer[yy+1:]

        time_recovery_samples = np.zeros(len(matched_filtered_samples), dtype=np.complex64)
        for i in range(0, len(matched_filtered_samples)):
            # push sample into interpolating filter
            self.tr_buf[:-1] = self.tr_buf[1:]
            self.tr_buf[-1] = matched_filtered_samples[i]

            # interpolate matched filter output
            self.tr_hi  = np.sinc(self.tr_th - self.tr_tau) * self.tr_w  # interpolating filter coefficients
            if i == 0:
                self.tr_hi /= np.sum(self.tr_hi)
            time_recovery_samples[i] = np.dot(self.tr_buf, self.tr_hi[::-1])    # compute matched filter output

            # take (approximate) derivative of filter output
            self.tr_buf_dz[:-1] = self.tr_buf_dz[1:]
            self.tr_buf_dz[-1] = time_recovery_samples[i]
            dz = -np.dot(self.tr_buf_dz, np.array([-1, 0, 1]))

            # determine if an output sample needs to be computed
            self.tr_counter += 1
            if self.tr_counter >= self.samples_per_symbol:
                # decrement counter by samples per symbol
                self.tr_counter -= self.samples_per_symbol

                # compute timing error signal, accounting for delay
                if i == 0:
                    self.err = np.tanh( (dz * np.conj(self.tr_last_sample)).real )
                else:
                    self.err = np.tanh( (dz * np.conj(time_recovery_samples[i-1])).real )

                # update timing rate change
                self.tr_dtau += self.tr_alpha * self.err
                self.tr_tau  += self.tr_beta * self.err

            # update timing error
            self.tr_tau += self.tr_dtau / self.samples_per_symbol

        self.tr_last_sample = time_recovery_samples[-1]

        phase_comp_samples = np.zeros(len(time_recovery_samples), dtype=np.complex64)

        # Costas loop
        for idx, sample in enumerate(time_recovery_samples):
            signal_out = sample * np.exp(-1.0j * self.cl_phase_est)
            phase_comp_samples[idx] = signal_out
            phase_error = np.sign(signal_out.real)*signal_out.imag - np.sign(signal_out.imag)*signal_out.real
            self.cl_frequency_out += self.cl_beta * phase_error
            self.cl_phase_est += self.cl_alpha * phase_error + self.cl_frequency_out

        if self.symbol_offset == -1:
            self.symbol_offset = 0
            sym_offset_0 = np.sum(np.abs(phase_comp_samples[::self.samples_per_symbol]))
            sym_offset_1 = np.sum(np.abs(phase_comp_samples[1::self.samples_per_symbol]))
            if sym_offset_1 < sym_offset_0:
                self.symbol_offset = 1
        demod_symbols = phase_comp_samples[self.symbol_offset::self.samples_per_symbol]
        self.symbols = demod_symbols

        bits = np.zeros(len(demod_symbols))
        angles = np.zeros(len(demod_symbols))

        for xx in range(0, len(demod_symbols)):
            if xx == 0:
                angle = np.angle(demod_symbols[xx], deg=True) - np.angle(self.last_symbol, deg=True)
            else:
                angle = np.angle(demod_symbols[xx], deg=True) - np.angle(demod_symbols[xx-1], deg=True)
            if angle > 180:
                angle -= 360
            if angle < -180:
                angle += 360
            angles[xx] = angle
            bit = 0
            if angle > 0: bit = 1
            bits[xx] = bit

        pos_angles = angles[np.where(angles > 0)]
        neg_angles = angles[np.where(angles < 0)]
        if len(pos_angles) > 0:
            self.ave_angles_above_zero = np.mean(pos_angles)
        if len(neg_angles) > 0:
            self.ave_angles_below_zero = np.mean(neg_angles)

        self.last_symbol = demod_symbols[-1]
        self.bit_string += ''.join([str(int(bit)) for bit in bits])

        num_of_possible_packets = len(self.bit_string) / self.size_of_packets
        if num_of_possible_packets > 20 and self.bit_offset == -1:
            # print("Calculate possible bit offset")
            # for all bit offsets (of the length of the packets)
            # calculate a score for most valid headers of that offset
            # this also checks a bit-reversed score (in case my endianness is reversed)
            scores = np.zeros(self.size_of_packets)
            revscores = np.zeros(self.size_of_packets)
            for xx in range(0, self.size_of_packets):
                for yy in range(xx, len(self.bit_string)-xx-8, self.size_of_packets):
                    if self.bit_string[yy:yy+8][::-1] in self.packet_headers:
                        scores[xx] += 1
                    if self.bit_string[yy:yy+8] in self.packet_headers:
                        revscores[xx] += 1

            self.reverse = False
            max_score = max(np.max(scores), np.max(revscores))
            if max_score > 15 and float(max_score)/num_of_possible_packets > 0.5:
                if np.max(scores) < np.max(revscores):
                    self.reverse = True
                if self.reverse:
                    self.bit_offset = np.argmax(revscores)
                else:
                    self.bit_offset = np.argmax(scores)
                self.bit_string = self.bit_string[self.bit_offset:]
                print("Bit stream offset: {}".format(self.bit_offset))

        self.packets = []
        if len(self.bit_string) > 4 * self.size_of_packets and self.bit_offset != -1:
            bits_to_remove = 0
            last_packet_epheris = False
            for xx in range(0, len(self.bit_string) - 2 * self.size_of_packets, self.size_of_packets):
                if last_packet_epheris == True:
                    last_packet_epheris = False
                    continue
                packet = ''
                if self.reverse:
                    header = '{:02X}'.format(int(self.bit_string[xx:xx+8], 2))
                else:
                    header = '{:02X}'.format(int(self.bit_string[xx:xx+8][::-1], 2))

                ephemeris_header = self.packet_dict['Ephemeris']['hex_header']
                packet_length = 12*8
                if header == ephemeris_header:
                    packet_length = 24*8
                    last_packet_epheris = True

                for yy in range(0, packet_length, 8):
                    if len(self.bit_string) < xx+yy+8:
                        packet = ''
                        break
                    if self.reverse:
                        packet += '{:02X}'.format(int(self.bit_string[xx+yy:xx+yy+8], 2))
                    else:
                        packet += '{:02X}'.format(int(self.bit_string[xx+yy:xx+yy+8][::-1], 2))
                if packet != '':
                    bits_to_remove += packet_length
                    self.packets.append(packet)
            # self.bit_string = self.bit_string[bits_to_remove:]

        # If 90% or more of the packets don't pass the checksum
        # leave the bits in the bit string and reset the bit offset to -1
        # on the next loop the code will recheck the bit offset to see if a
        # better bit offset will allow us to decode the bits
        if len(self.packets) > 0:
            total_packets = len(self.packets)
            bad_packets = 0.0
            for packet in self.packets:
                if fletcher_checksum(packet) != '0000':
                    bad_packets += 1
            if float(bad_packets)/total_packets > 0.9:
                self.packets = []
                self.bit_offset = -1
                print("Packets failing checksum. Resetting bit offset.")
            else:
                self.bit_string = self.bit_string[bits_to_remove:]

        if num_of_possible_packets > 60 and self.bit_offset == -1:
            print("RESETTING!")
            self.init_default_values()

        packet_dict_list = []
        for packet in self.packets:

            # Compute the fletcher16 checksum over the whole packet
            # 0000 output is a good packet
            if fletcher_checksum(packet) != '0000':
                self.bad_packets += 1
                continue
            else:
                self.good_packets += 1
            if self.good_packets + self.bad_packets > 500:
                self.good_packets -= 0.5
                self.bad_packets -= 0.5
                if self.good_packets < 0: self.good_packets = 0.0
                if self.bad_packets < 0: self.bad_packets = 0.0

            packet_dict = {}
            for packet_type in self.packet_dict:
                packet_info = self.packet_dict[packet_type]
                if packet[:2] == packet_info['hex_header']:
                    packet_dict['packet_type'] = packet_type

                    for part, (start, stop) in packet_info['message_parts']:
                        packet_dict[part] = packet[start:stop]

                    if packet_type == 'Ephemeris':
                        payload = ''.join([packet[xx:xx+2] for xx in range(42, 2, -2)])

                        # calculate current satellite time
                        start_date = datetime(year=1980, month=1, day=6, hour=0, minute=0)
                        week_number = payload[:4]
                        time_of_week = payload[4:10]
                        this_week = start_date + timedelta(weeks=int(week_number, 16))
                        this_time = this_week + timedelta(seconds=int(time_of_week, 16))

                        # calculate satellite ECEF position
                        zdot = payload[10:15][::-1]
                        ydot = payload[15:20][::-1]
                        xdot = payload[20:25][::-1]
                        zpos = payload[25:30][::-1]
                        ypos = payload[30:35][::-1]
                        xpos = payload[35:40][::-1]

                        max_r_sat = 8378155.0
                        val_20_bits = 1048576.0

                        x_temp = int(xpos[:2][::-1], 16) + 256. * int(xpos[2:4][::-1], 16) + 256**2 * int(xpos[4:], 16)
                        x_ecef = ((2*x_temp*max_r_sat)/val_20_bits - max_r_sat)
                        y_temp = int(ypos[:2][::-1], 16) + 256. * int(ypos[2:4][::-1], 16) + 256**2 * int(ypos[4:], 16)
                        y_ecef = ((2*y_temp*max_r_sat)/val_20_bits - max_r_sat)
                        z_temp = int(zpos[:2][::-1], 16) + 256. * int(zpos[2:4][::-1], 16) + 256**2 * int(zpos[4:], 16)
                        z_ecef = ((2*z_temp*max_r_sat)/val_20_bits - max_r_sat)

                        lat, lon, alt = ecef_to_lla(x_ecef, y_ecef, z_ecef)
                        self.sat_lat = lat
                        self.sat_lon = lon
                        self.sat_alt = alt
                        packet_dict['lat'] = lat
                        packet_dict['lon'] = lon
                        packet_dict['alt'] = alt
                        packet_dict['week_num'] = int(week_number, 16)
                        packet_dict['gps_time'] = this_time
                    break

            # Unrecognized just means I don't know what these packets are for
            # would also happen if the header is corrupted
            if 'packet_type' not in packet_dict:
                packet_dict['packet_type'] = 'Unrecognized'
                packet_dict['data'] = packet[2:]
        return packet_dict_list

    def parse_packets(self):
        for packet in self.packets:
            output = ''

            # Compute the fletcher16 checksum over the whole packet
            # 0000 output is a good packet
            if fletcher_checksum(packet) != '0000':
                self.bad_packets += 1
                output += '### '
                continue
            else:
                self.good_packets += 1
            if self.good_packets + self.bad_packets > 500:
                self.good_packets -= 0.5
                self.bad_packets -= 0.5
                if self.good_packets < 0: self.good_packets = 0.0
                if self.bad_packets < 0: self.bad_packets = 0.0
            for packet_type in self.packet_dict:
                packet_info = self.packet_dict[packet_type]
                if packet[:2] == packet_info['hex_header']:
                    output += '{}: '.format(packet_type)
                    for part, (start, stop) in packet_info['message_parts']:
                        output += '{}: {} '.format(part, packet[start:stop])
                    print(output)
                    if packet_type in ['Sync', 'Ephemeris']:
                        for part, (start, stop) in packet_info['message_parts']:
                            if part == 'Sat ID':
                                self.sat_id = packet[start:stop]
                    if packet_type == 'Ephemeris':
                        payload = ''.join([packet[xx:xx+2] for xx in range(42, 2, -2)])

                        # calculate current satellite time
                        start_date = datetime(year=1980, month=1, day=6, hour=0, minute=0)
                        week_number = payload[:4]
                        time_of_week = payload[4:10]
                        this_week = start_date + timedelta(weeks=int(week_number, 16))
                        this_time = this_week + timedelta(seconds=int(time_of_week, 16))
                        print("\tCurrent satellite time: {} Z".format(this_time))

                        # calculate satellite ECEF position
                        zdot = payload[10:15][::-1]
                        ydot = payload[15:20][::-1]
                        xdot = payload[20:25][::-1]
                        zpos = payload[25:30][::-1]
                        ypos = payload[30:35][::-1]
                        xpos = payload[35:40][::-1]

                        max_r_sat = 8378155.0
                        val_20_bits = 1048576.0

                        x_temp = int(xpos[:2][::-1], 16) + 256. * int(xpos[2:4][::-1], 16) + 256**2 * int(xpos[4:], 16)
                        x_ecef = ((2*x_temp*max_r_sat)/val_20_bits - max_r_sat)
                        y_temp = int(ypos[:2][::-1], 16) + 256. * int(ypos[2:4][::-1], 16) + 256**2 * int(ypos[4:], 16)
                        y_ecef = ((2*y_temp*max_r_sat)/val_20_bits - max_r_sat)
                        z_temp = int(zpos[:2][::-1], 16) + 256. * int(zpos[2:4][::-1], 16) + 256**2 * int(zpos[4:], 16)
                        z_ecef = ((2*z_temp*max_r_sat)/val_20_bits - max_r_sat)

                        lat, lon, alt = ecef_to_lla(x_ecef, y_ecef, z_ecef)
                        self.sat_lat = lat
                        self.sat_lon = lon
                        self.sat_alt = alt
                        print("\tLat/Lon:       {:8.4f}, {:8.4f}, Altitude: {:6.1f} km".format(lat, lon, alt/1000.0))
                    break

            # Unrecognized just means I don't know what these packets are for
            # would also happen if the header is corrupted
            if output in ['', '### ']:
                print("{}Unrecognized packet: {}".format(output, packet))