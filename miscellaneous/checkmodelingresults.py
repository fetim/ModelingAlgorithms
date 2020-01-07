#!/usr/bin/python

import os
import numpy as np

f = open('../parameters/2D_acoustic_modeling.dat', 'r') # 'r' = read
parameters = np.genfromtxt(f,delimiter='')

Nx     = int(parameters[0])
Nz     = int(parameters[1])
dx     = parameters[2]
dz     = parameters[3]
Nt     = int(parameters[4])
dt     = parameters[5]
Nshot  = parameters[6]
fcut   = parameters[7]

print("Nx    = ",Nx  )
print("Nz    = ",Nz  )
print("dx    = ",dx  )
print("dz    = ",dz  )
print("Nt    = ",Nt  )
print("dt    = ",dt  )
print("Nshot = ",Nshot)
print("fcut  = ",fcut)


folder = "../fortran/openacc/"
# Check seismogram with seismic unix
command = " ximage n1=" + str(Nt) +" < " + folder + "seismogram.bin &"
print(command)
os.system(command)

# Check snapshot with seismic unix
command = " ximage n1=" + str(Nz) + " < " + folder + "snapshots.bin &"
# command = " xmovie n1=" + str(Nz) + " n2=" + str(Nx) + " sleep=1 loop=2 < " + folder + "snapshots.bin &"

print(command)
os.system(command)