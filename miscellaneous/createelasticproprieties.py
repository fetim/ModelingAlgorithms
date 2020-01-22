import numpy as np
import matplotlib.pyplot as pl


# parameters
name = "modelsBuzios/xline_3825_buzios_1001nz_2401nx_dz10m_dx25m.bin"
Nx = 2401
Nz = 1001
dx = 25
dz = 10

trace = np.arange(0,dx*Nx,dx)    
depth = np.arange(0,dz*Nz,dz)    

# import velocity model
velocity_horizon = qc.readbinaryfile(Nz,Nx,name)

##### define Vs
velocity_horizon_vs = velocity_horizon*0.5

# Edit water layer for vs
print("filling layer", layer)
for ii in range(0,Nx):
    for jj in range(0,Nz):
        # first layer
        layer = 1 
        if (jj <= horizons[ii,layer]):                    
            velocity_horizon_vs[jj,ii] = 0

# plot in same window
fig, ax = pl.subplots()
im = ax.imshow(velocity_horizon_vs,cmap='jet')

for layer in range(1,N_layers):                
        ax.plot(horizons[:,0],horizons[:,layer])     

fig.colorbar(im,ax=ax)
pl.show(block=False) 

# define rho (Gardner)                
velocity_horizon_rho = 310*np.power(velocity_horizon,0.25)

# Edit water layer for rho
print("filling layer", layer)
for ii in range(0,Nx):
    for jj in range(0,Nz):
        # first layer
        layer = 1 
        if (jj <= horizons[ii,layer]):                    
            velocity_horizon_rho[jj,ii] = 1000

fig, ax = pl.subplots()
im = ax.imshow(velocity_horizon_rho,cmap='jet')

for layer in range(1,N_layers):            
        ax.plot(horizons[:,0],horizons[:,layer])     

fig.colorbar(im,ax=ax)
pl.show(block=False) 

# # write velocity S
# outfile = name + '_layer_cake_vs.bin'
# qc.savebinaryfile(Nz,Nx,velocity_horizon_vs,outfile)

# # write density
# outfile = name + '_layer_cake_rho.bin'
# qc.savebinaryfile(Nz,Nx,velocity_horizon_rho,outfile)
