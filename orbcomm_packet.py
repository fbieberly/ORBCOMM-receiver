# List of the orbcomm satellites that I currently know are transmitting



packet_headers = [
                '01100101',     # 65 A8 F9, Synchronization
                '00011010',     # 1A 30/31/32, Single/Multi-line message packet
                '00011011',     # 1B 10, Uplink information in 1 packet format
                '00011100',     # 1C 30/31/32, Downlink information in 3 packet format
                '00011101',     # 1D 10, Network Control Centre information in 1 packet format
                '00011110',     # 1E, Fill packet
                '00011111',     # 1F, Ephemeris in 2 packet format
                '00100010',     # 22, Orbital elements in 1 packet format
                    ]

new_packet_headers = [
                        '00000011', # 03, 
                        '00111100', # 3C, fill packet
                        '10111100', # BC, fill packet? (Lots of repeated packets)
                        '01010010',
                        '01011000', # 58 0F (0101100000001111)
                        '10110000', # D8 0F (1101100000001111)
                        '01100010',
                        '01101111',
                        '10001000',
                        '10100000',
                        '11000000',
                        '11011110',
                        '11110001',
                        ]