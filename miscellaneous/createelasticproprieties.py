import numpy as np
import matplotlib.pyplot as pl

import qualitycontrolfunctions as qc

# parameters
name = "modelsBuzios/xline_3825_buzios_1001nz_2401nx_dz10m_dx25m.bin_layer_cake_from_10horizons.bin"
Nx = 2401
Nz = 1001
dx = 25
dz = 10

trace = np.arange(0,dx*Nx,dx)    
depth = np.arange(0,dz*Nz,dz)    

# load horizons
fundomar = np.loadtxt("2D_horizons/Fundo_do_Mar_XL3825_2D_davidff.txt")
pos_sal1 = np.loadtxt("2D_horizons/Pos-Sal_1_XL3825_2D_davidff.txt")
pos_sal2 = np.loadtxt("2D_horizons/Pos-Sal_2_XL3825_2D_davidff.txt")
pos_sal3 = np.loadtxt("2D_horizons/Pos-Sal_3_XL3825_2D_davidff.txt")
toposal  = np.loadtxt("2D_horizons/Topo_Sal_XL3825_2D_davidff.txt")
basesal  = np.loadtxt("2D_horizons/Base_Sal_XL3825_2D_davidff.txt")

N_horizon = 7
N_layers = N_horizon + 1
horizons = np.zeros([Nx,N_layers])

# Create horizon table
horizons[:,0]  = np.arange(0,Nx) # trace index
horizons[:,1]  = np.rint(fundomar[:-2,4]/dz -1)
horizons[:,2]  = np.rint(pos_sal1[:-1,4]/dz -1)
horizons[:,3]  = np.rint(pos_sal2[:,4]/dz -1)
horizons[:,4]  = np.rint(pos_sal3[:-1,4]/dz -1)
horizons[:,5]  = np.rint(toposal[:,4]/dz -1)
horizons[:,6]  = np.rint(basesal[:,4]/dz -1)

# repeat last horizon for base of reservoir
horizons[:,7]  = np.rint(basesal[:,4]/dz -1) + 55

# import velocity model
velocity_P = qc.readbinaryfile(Nz,Nx,name)

# plot in same window
fig, ax = pl.subplots()
im = ax.imshow(velocity_P,cmap='jet')

for layer in range(1,N_layers):
        # pl.plot(horizon[:,0,layer],horizon[:,1,layer])
        ax.plot(horizons[:,0],horizons[:,layer])     

fig.colorbar(im,ax=ax)
pl.title("Input P velocity model")
pl.show(block=False) 

# Vp - Vs Relation
# Castanha relation - Rosa (2010) apud Castanha et al. pag.513    vs = ai + bivp + civp^2

a     = np.zeros([N_layers,1])
a[0]  = -1172
a[1]  = -1172
a[2]  = -1172
a[3]  = -1172
a[4]  = -1172
a[5]  = -1172
a[6]  = -1172
a[7]  = -1172

b     = np.zeros([N_layers,1])
b[0]  = 0.862
b[1]  = 0.862
b[2]  = 0.862
b[3]  = 0.862
b[4]  = 0.862
b[5]  = 0.862
b[6]  = 0.862
b[7]  = 0.862

velocity_S = np.zeros([Nz,Nx])

# first layer
layer = 1 
print("Creating VS model")
print("filling layer", layer)
for ii in range(0,Nx):
    for jj in range(0,Nz):        
        if (jj <= horizons[ii,layer]):                    
            velocity_S[jj,ii] = 0

# 2nd layer - (N-1)th layer
for layer in range(2,N_layers):
    print("filling layer", layer)
    for ii in range(0,Nx):
        for jj in range(0,Nz):   
            #layer =2
            if (jj > horizons[ii,layer-1] and jj <= horizons[ii,layer]):                    
                velocity_S[jj,ii] = velocity_P[jj,ii]*b[layer-1] + a[layer-1]

# last layer
layer = N_layers
print("filling layer", layer)
for ii in range(0,Nx):
    for jj in range(200,Nz):
        if (jj > horizons[ii,layer-1]):
            velocity_S[jj,ii] = velocity_P[jj,ii]*b[layer-1] + a[layer-1]

# plot in same window
fig, ax = pl.subplots()
im = ax.imshow(velocity_S,cmap='jet')

# for layer in range(1,N_layers):                
#         ax.plot(horizons[:,0],horizons[:,layer])     
pl.title("Shear Velocity model")
fig.colorbar(im,ax=ax)
pl.show(block=False) 

# VP - RHO relations
# Gardner relation - Rosa (2010) apud Gardner et al. (1974)  pag. 496rho = a*vp^{b}
# 1 feet = 0.3048 m
# convert g/cm^3 > kg/m^3 and feet/s => m/s

a     = np.zeros([N_layers,1])
a[0]  = 0.23
a[1]  = 0.23
a[2]  = 0.23
a[3]  = 0.23
a[4]  = 0.23
a[5]  = 0.23
a[6]  = 0.23
a[7]  = 0.23

b     = np.zeros([N_layers,1])
b[0]  = 0.25
b[1]  = 0.25
b[2]  = 0.25
b[3]  = 0.25
b[4]  = 0.25
b[5]  = 0.25
b[6]  = 0.25
b[7]  = 0.25

rho = np.zeros([Nz,Nx])
# first layer
layer = 1 
print("Creating density model")
print("filling layer", layer)
for ii in range(0,Nx):
    for jj in range(0,Nz):        
        if (jj <= horizons[ii,layer]):                    
            rho[jj,ii] = 1000

# 2nd layer - (N-1)th layer
for layer in range(2,N_layers):
    print("filling layer", layer)
    for ii in range(0,Nx):
        for jj in range(0,Nz):   
            #layer =2
            if (jj > horizons[ii,layer-1] and jj <= horizons[ii,layer]):                    
                rho[jj,ii] = (a[layer-1]*np.power(velocity_P[jj,ii]/0.3048,b[layer-1]))*1000

# last layer
layer = N_layers
print("filling layer", layer)
pl.title("Density model")
for ii in range(0,Nx):
    for jj in range(200,Nz):
        if (jj > horizons[ii,layer-1]):
            rho[jj,ii] = (a[layer-1]*np.power(velocity_P[jj,ii]/0.3048,b[layer-1]))*1000
            
# plot in same window
fig, ax = pl.subplots()
im = ax.imshow(rho,cmap='jet')

# for layer in range(1,N_layers):                
#         ax.plot(horizons[:,0],horizons[:,layer])     

fig.colorbar(im,ax=ax)
pl.show(block=False) 

# write velocity S
outfile = name + '_vs.bin'
qc.savebinaryfile(Nz,Nx,velocity_S,outfile)

# write density
outfile = name + '_rho.bin'
qc.savebinaryfile(Nz,Nx,rho,outfile)
