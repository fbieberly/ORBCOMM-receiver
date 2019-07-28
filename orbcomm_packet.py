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
                                            ('Code', (0, 6)),  # Part name, (start idx, stop idx)
                                            ('Sat ID', (6, 8)),  
                                            ],
                            },
                    'Message':{
                            'header':'00011010',
                            'hex_header':'1A',
                            'message_parts':[
                                            ('Part', (3, 4)),
                                            ('Total length', (2, 3)),
                                            ('Data', (4, 22)),
                                            ],
                            },
                    'Uplink_info':{
                            'header':'00011011',
                            'hex_header':'1B',
                            'message_parts':[
                                            ('Total length', (2, 3)),
                                            ('Part', (3, 4)),
                                            ('Data', (4, 22)),
                                            ],

                            },
                    'Downlink_info':{
                            'header':'00011100',
                            'hex_header':'1C',
                            'message_parts':[
                                            ('Total length', (2, 3)),
                                            ('Part', (3, 4)),
                                            ('Data', (4, 22)),
                                            ],

                            },
                    'Network':{
                            'header':'00011101',
                            'hex_header':'1D',
                            'message_parts':[
                                            ('Total length', (2, 3)),
                                            ('Part', (3, 4)),
                                            ('Data', (4, 22)),
                                            ],

                            },
                    'Fill':{
                            'header':'00011110',
                            'hex_header':'1E',
                            'message_parts':[ ('Data', (2, 22)),
                                            ],

                            },
                    'Ephemeris':{
                            'header':'00011111',
                            'hex_header':'1F',
                            'message_parts':[
                                            ('Sat ID', (2, 4)),  
                                            ('Data', (4, 48)),  
                                            ],

                            },
                    'Oribital':{
                            'header':'00100010',
                            'hex_header':'22',
                            'message_parts':[
                                            ],

                            },
                    }