#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computes total hours since 1/1/1900 00:00:00 gregorian calendar
"""

from datetime import datetime, timedelta

def total_hours(day, month, year, hour):
    di = datetime(year= 1900, month = 1, day = 1, hour = 0)
    df = datetime(year= year, month = month, day = day, hour = hour)
    return int((df - di).total_seconds()/3600)
