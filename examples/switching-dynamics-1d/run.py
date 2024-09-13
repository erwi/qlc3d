#!/bin/python

import os
import fnmatch
import shutil
import subprocess
import math
import matplotlib.pyplot as plt

executable = "../../cmake-build-debug/qlc3d/qlc3d"


# delete old results, if any exist
shutil.rmtree("res", ignore_errors=True)

# run the simulation
subprocess.call([executable, "settings.txt"])
    
# read the results from the files output in the "res directory" 
os.chdir('res')
result_files = fnmatch.filter(os.listdir(), "dirstacksz*.csv")
result_files.sort()
print("result files = " + str(result_files));


max_tilts = []
times = []
counter = 0;
save_time = 1e-4; # must be same valueas SaveTime in settings file
for res in result_files:
    fid = open(res)
    line = fid.readline() # skip first line
    line = fid.readline() # data in second line
    fid.close();
    
    # convert the line to arrays with director x,y,z components
    splits = [float(v) for v in line.split(',')]
    n_x = splits[0:30:3]
    n_y = splits[1:30:3]
    n_z = splits[2:30:3]
    # convert directo to tilt angle
    tilts = [math.degrees(math.asin(math.fabs(nz))) for nz in n_z]
    
    # mid-plane tilt angle for current time step
    max_tilts.append(max(tilts))
    times.append(counter * save_time)
    counter += 1;

# plot the tilt angle as function of time.
plt.plot(times, max_tilts, '.-')
plt.xlabel("time [seconds]")
plt.ylabel("max tilt angle [degrees]")
plt.show()