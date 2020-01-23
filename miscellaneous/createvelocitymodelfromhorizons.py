import numpy as np
import matplotlib.pyplot as pl
import qualitycontrolfunctions as qc

# parameters
Nx = 2401
Nz = 1001
dx = 25
dz = 10

trace = np.arange(0,dx*Nx,dx)    
depth = np.arange(0,dz*Nz,dz)   

# reference velocity model
name = "modelsBuzios/xline_3825_buzios_1001nz_2401nx_dz10m_dx25m.bin"

# load horizons
fundomar = np.loadtxt("2D_horizons/Fundo_do_Mar_XL3825_2D_davidff.txt")
pos_sal1 = np.loadtxt("2D_horizons/Pos-Sal_1_XL3825_2D_davidff.txt")
pos_sal2 = np.loadtxt("2D_horizons/Pos-Sal_2_XL3825_2D_davidff.txt")
pos_sal3 = np.loadtxt("2D_horizons/Pos-Sal_3_XL3825_2D_davidff.txt")
toposal  = np.loadtxt("2D_horizons/Topo_Sal_XL3825_2D_davidff.txt")
basesal  = np.loadtxt("2D_horizons/Base_Sal_XL3825_2D_davidff.txt")

N_horizon = 11
N_layers = N_horizon + 1
horizons = np.zeros([Nx,N_layers])

# Create horizon table
horizons[:,0]  = np.arange(0,Nx) # trace index
horizons[:,1]  = np.rint(fundomar[:-2,4]/dz -1)

horizons[:,3]  = np.rint(pos_sal1[:-1,4]/dz -1)
# creating horizon between fundomar and pos_sal1
horizons[:,2] = np.rint((horizons[:,1] + horizons[:,3])/2)

horizons[:,4]  = np.rint(pos_sal2[:,4]/dz -1)

horizons[:,6]  = np.rint(pos_sal3[:-1,4]/dz -1)
# creating horizon between pos_sal2 and pos_sal3
horizons[:,5] = np.rint((horizons[:,4] + horizons[:,6])/2)

horizons[:,9]  = np.rint(toposal[:,4]/dz -1)

# creating horizon between pos_sal3 and toposal
horizons[:,7] = np.rint((horizons[:,6] + horizons[:,9])/2)

# creating horizon between pos_sal3.1 and toposal
horizons[:,8] = np.rint((horizons[:,7] + horizons[:,9])/2)

horizons[:,10]  = np.rint(basesal[:,4]/dz -1)

# repeat last horizon for base of reservoir
horizons[:,11]  = np.rint(basesal[:,4]/dz -1) + 55


# set velocity for each layer        
horizon_vel     = np.zeros([N_layers,1])
horizon_vel[0]  = 1500
horizon_vel[1]  = 2050
horizon_vel[2]  = 2750
horizon_vel[3]  = 3000
horizon_vel[4]  = 3050
horizon_vel[5]  = 3150
horizon_vel[6]  = 3200
horizon_vel[7]  = 3500
horizon_vel[8]  = 3800
horizon_vel[9]  = 4500
horizon_vel[10] = 3800
horizon_vel[11] = 4600



velocity_horizon = np.zeros([Nz,Nx])

# first layer
layer = 1 
print("filling layer", layer)
for ii in range(0,Nx):
    for jj in range(0,Nz):        
        if (jj <= horizons[ii,layer]):                    
            velocity_horizon[jj,ii] = horizon_vel[layer-1]                

# 2nd layer - (N-1)th layer
for layer in range(2,N_layers):
    print("filling layer", layer)
    for ii in range(0,Nx):
        for jj in range(0,Nz):   
            #layer =2
            if (jj > horizons[ii,layer-1] and jj <= horizons[ii,layer]):                    
                velocity_horizon[jj,ii] = horizon_vel[layer-1]  

# last layer
layer = N_layers
print("filling layer", layer)
for ii in range(0,Nx):
    for jj in range(0,Nz):
        if (jj > horizons[ii,layer-1]):
            velocity_horizon[jj,ii] = horizon_vel[layer-1]

# plot in same window
fig, ax = pl.subplots()
im = ax.imshow(velocity_horizon,cmap='jet')

for layer in range(1,N_layers):
        # pl.plot(horizon[:,0,layer],horizon[:,1,layer])
        ax.plot(horizons[:,0],horizons[:,layer])     

fig.colorbar(im,ax=ax)
pl.show(block=False) 

# write velocity P
outfile = name + '_layer_cake_from_horizons.bin'
qc.savebinaryfile(Nz,Nx,velocity_horizon,outfile)