#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert from time from WRF to time from ECMWF
"""

import datetime
import numpy as np
from re import split

def convert_to_hours_since_1990(time, time_units):
    
    ref_time = split('since ', time_units)
    ref_time = datetime.datetime.strptime(ref_time[1], '%Y-%m-%d %H:%M:%S')
    diff_hours = ref_time - datetime.datetime(1900,1,1,0)
    diff_hours = diff_hours.total_seconds()/3600
    
    times = time/60
    times = times + diff_hours
    return times
    
    
    

