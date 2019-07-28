##############################################################################
#
# Author: Frank Bieberly
# Date: 17 July 2019
# Name: sat_db.py
# Description: 
# This is a list of the orbcomm satellites that I'm pretty sure are still
# transmitting.
# 
#
##############################################################################


active_orbcomm_satellites = {
            # No signal in several passes:
                # FM08, FM09, FM06, FM07, FM04, FM05, FM02, FM03, FM01, FM19,
                # FM18, FM11, FM10, FM12, FM13, FM15, FM14, FM17, FM16, FM20,
                # FM21, FM22, FM23, FM24, FM25, FM26, FM28, FM32, FM31, FM30, 
                # FM35, 

            # Not enough info to know: 
                # FM104, FM107, FM106, FM105, FM27, FM29, FM33, FM37, FM36, FM34,
                # FM111, FM119, FM115
            'orbcomm fm103':{
                            'norad_id':40091,
                            'frequencies':[137250000.0, 137312500.0]
                            },
            'orbcomm fm107':{
                            'norad_id':40087,
                            'frequencies':[137250000.0, 137312500.0]
                            },
            'orbcomm fm108':{
                            'norad_id':41187,
                            'frequencies':[137460000.0, 137712500.0]
                            },
            'orbcomm fm109':{
                            'norad_id':40086,
                            'frequencies':[137250000.0, 137312500.0]
                            },
            'orbcomm fm110':{
                            'norad_id':41182,
                            'frequencies':[137287500.0, 137737500.0]
                            },
            'orbcomm fm112':{
                            'norad_id':41184,
                            'frequencies':[137662500.0, 137800000.0]
                            },
            'orbcomm fm113':{
                            'norad_id':41185,
                            'frequencies':[137662500.0, 137800000.0]
                            },
            'orbcomm fm114':{
                            'norad_id':41179,
                            'frequencies':[137287500.0, 137737500.0]
                            },
            'orbcomm fm116':{
                            'norad_id':41189,
                            'frequencies':[137662500.0, 137800000.0]
                            },
            'orbcomm fm117':{
                            'norad_id':41188,
                            'frequencies':[137460000.0, 137712500.0]
                            },
            'orbcomm fm118':{
                            'norad_id':41183,
                            'frequencies':[137287500.0, 137737500.0]
                            }
            }