import os
from calendar import monthrange
import numpy as np

month1 = 3
month2 = 7
mms = np.arange(month1, month2+1)


for month in mms:
    dds = monthrange(2015, month)[1]
    for day in range(1, dds+1):
        os.system('mv wrf_2015-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_00:00:00 wrfout_d01_2015-'+str(month).zfill(2)+'-'+str(day).zfill(2)+'_00:00:00')
