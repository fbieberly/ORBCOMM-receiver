##############################################################################
#
# Author: Frank Bieberly
# Date: 16 July 2019
# Name: orbcomm_packet.py
# Description:
# Dictionary of information about orbcomm packets.
#
#
##############################################################################

packet_dict = {
                    'Sync':{
                            'header':'01100101',
                            'hex_header':'65',
                            'message_parts':[
                                            ('code', (0, 6)),  # Part name, (start idx, stop idx)
                                            ('sat_id', (6, 8)),
                                            ],
                            },
                    'Message':{
                            'header':'00011010',
                            'hex_header':'1A',
                            'message_parts':[
                                            ('msg_packet_num', (3, 4)),
                                            ('msg_total_length', (2, 3)),
                                            ('data', (4, 22)),
                                            ],
                            },
                    'Uplink_info':{
                            'header':'00011011',
                            'hex_header':'1B',
                            'message_parts':[
                                            ('msg_packet_num', (3, 4)),
                                            ('msg_total_length', (2, 3)),
                                            ('data', (4, 22)),
                                            ],

                            },
                    'Downlink_info':{
                            'header':'00011100',
                            'hex_header':'1C',
                            'message_parts':[
                                            ('msg_packet_num', (3, 4)),
                                            ('msg_total_length', (2, 3)),
                                            ('data', (4, 22)),
                                            ],

                            },
                    'Network':{
                            'header':'00011101',
                            'hex_header':'1D',
                            'message_parts':[
                                            ('msg_packet_num', (3, 4)),
                                            ('msg_total_length', (2, 3)),
                                            ('data', (4, 22)),
                                            ],

                            },
                    'Fill':{
                            'header':'00011110',
                            'hex_header':'1E',
                            'message_parts':[ ('data', (2, 22)),
                                            ],

                            },
                    'Ephemeris':{
                            'header':'00011111',
                            'hex_header':'1F',
                            'message_parts':[
                                            ('sat_id', (2, 4)),
                                            ('data', (4, 48)),
                                            ],

                            },
                    'Oribital':{
                            'header':'00100010',
                            'hex_header':'22',
                            'message_parts':[
                                            ],

                            },

                    # Currently unknown packet types
                    'Unknown 1':{
                            'header':'00001010',
                            'hex_header':'0A',
                            'message_parts':[
                                            ('header', (0, 2)),
                                            ('data', (2, 48)),
                                            ],

                            },
                    'Unknown 2':{
                            'header':'00001011',
                            'hex_header':'0B',
                            'message_parts':[
                                            ('header', (0, 2)),
                                            ('data', (2, 48)),
                                            ],

                            },
                    'Unknown 3':{
                            'header':'00001101',
                            'hex_header':'0D',
                            'message_parts':[
                                            ('header', (0, 2)),
                                            ('data', (2, 48)),
                                            ],

                            },
                    'Unknown 4':{
                            'header':'00001110',
                            'hex_header':'0E',
                            'message_parts':[
                                            ('header', (0, 2)),
                                            ('data', (2, 48)),
                                            ],

                            },
                    'Unknown 5':{
                            'header':'00010011',
                            'hex_header':'13',
                            'message_parts':[
                                            ('header', (0, 2)),
                                            ('data', (2, 48)),
                                            ],

                            },
                    'Unknown 6':{
                            'header':'00010100',
                            'hex_header':'18',
                            'message_parts':[
                                            ('header', (0, 2)),
                                            ('data', (2, 48)),
                                            ],

                            },
                    }
