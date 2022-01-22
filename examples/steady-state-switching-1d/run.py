#!/bin/python

import os
import subprocess
import math
import matplotlib.pyplot as plt

executable = "../../cmake-build-release/qlc3d/qlc3d"

potentials = [p / 10 for p in range(31)]
tilts = []

referencePotentials = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3]
referenceTilts = [5, 5.71, 9.2, 24.0, 44.8, 60.0, 70.46, 81.8, 86.6, 88.6]


for pot in potentials:
    print("potential=" + str(pot))
    os.environ["QLC3D_POTENTIAL_1"] = str(pot)
    subprocess.call([executable, "settings.txt"])
    
    # read the director from the last result file in the results directory
    fid = open("res/dirstacksz-final.csv")
    line = fid.readline()
    line = fid.readline()
    fid.close();
    

    splits = line.split(',')
    n_z = float(splits[2])
    tilt = math.degrees(math.asin(n_z))

    print("mid plane director " + str(line) + " tilt [degrees] = " + str(tilt))
    tilts.append(tilt)

print("tilt angles = " + str(tilts))
plt.plot(potentials, tilts, label = "qlc3d")
plt.plot(referencePotentials, referenceTilts, ".", label = "lc3k")
plt.legend()
plt.xlabel("potential [V.]")
plt.ylabel("mid-plane tilt [degrees]")
plt.show()