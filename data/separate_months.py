import os 
from calendar import monthrange





files = os.listdir('./MET/')
orig_file = files[0]

os.chdir('./MET')

first_time = 0
for i in range(12):
    days = monthrange(2015, i+1)[1]
    
    last_time = first_time + days*24 -1
    out_file = 'ERA5_'+str(i+1).zfill(2) +'_2015.nc'
    os.system('ncks -d time,'+str(first_time)+','+str(last_time) + ' '+orig_file + ' '+ out_file)
    first_time += days*24
