#!/opt/anaconda3/bin/python
'''
This scricpt create a 2D velocity model from given horizons
The velocities for each layer can be defined by user or
based on mean of another (e.g. smooth) velocity model.
'''
import numpy as np
import matplotlib.pyplot as pl
import qualitycontrolfunctions as qc
from mpl_toolkits.axes_grid1 import make_axes_locatable # colobar adjust

# parameters
folder ="/mnt/Desktop/Dados_Buzios/" 
name = "xline_3824_buzios_1001z_2403x_dz10m_dx25m_interp"
outfile = folder + name +'_Lcake_horiz12_mean_vag_velmig'

Nx = 2403
Nz = 1001
dx = 25
dz = 10

trace = np.arange(0,dx*Nx,dx)    
depth = np.arange(0,dz*Nz,dz)    


# import velocity model
velocity = qc.readbinaryfile(Nz,Nx,folder+name+".bin")

# load horizons
folder_horizon = folder+"2D_Horizons_IL4328_XL3824/"
fundomar = np.loadtxt(folder_horizon+"Fundo_do_Mar_XL3824_2D_davidff.txt")
pos_sal1 = np.loadtxt(folder_horizon+"Pos-sal_1_XL3824_2D_davidff.txt")
pos_sal2 = np.loadtxt(folder_horizon+"Pos-sal_2_XL3824_2D_raissamss.txt")
pos_sal3 = np.loadtxt(folder_horizon+"Pos-sal_3_XL3824_2D_davidff.txt")
toposal  = np.loadtxt(folder_horizon+"Topo_Sal_XL3824_2D_davidff.txt")
basesal  = np.loadtxt(folder_horizon+"Base_Sal_XL3824_2D_davidff.txt")


# number=1 ## prepare horizons and create one horizons in reservoir
number=2## create fake horizons between true horizons    
if (number == 1): ## prepare horizons and create one horizons in reservoir
    N_horizon = 8
    N_layers = N_horizon + 1
    horizons = np.zeros([Nx,N_layers])    
    ### Create horizon table    
    # trace index
    horizons[:,0]  = np.arange(0,Nx) 

    horizons[:,1]  = np.rint(fundomar[:,4]/dz -1)
    horizons[:,2]  = np.rint(pos_sal1[:,4]/dz -1)
    horizons[:,3]  = np.rint(pos_sal2[:,4]/dz -1)
    horizons[:,4]  = np.rint(pos_sal3[:,4]/dz -1)
    horizons[:,5]  = np.rint(toposal[:,4]/dz -1)
    horizons[:,6]  = np.rint(basesal[:,4]/dz -1)
    # repeat base of salt horizon for layer inside reservoir salt's base + 170m
    horizons[:,7]  = np.rint(qc.smooth1D(basesal[:,4],50)/dz -1) + 17
    # repeat base of salt horizon for base of reservoir at 370 m bellow 
    horizons[:,8]  = np.rint(qc.smooth1D(basesal[:,4],100)/dz -1) + 37

    ## set velocity for each layer
    # mean velocity from 2anp borehole        
    horizon_vel     = np.zeros([N_layers,1])
    horizon_vel[0]  = 1500
    horizon_vel[1]  = 1800
    horizon_vel[2]  = 1950
    horizon_vel[3]  = 2200
    horizon_vel[4]  = 2600
    horizon_vel[5]  = 4500
    horizon_vel[6]  = 5450
    horizon_vel[7]  = 4500
    horizon_vel[8]  = 4800
    
elif (number ==2): ## create fake horizons between true horizons    

    N_horizon = 12
    N_layers = N_horizon + 1
    horizons = np.zeros([Nx,N_layers])
    ### Create horizon table    
    # trace index
    horizons[:,0]  = np.arange(0,Nx)

    horizons[:,1]  = np.rint(fundomar[:,4]/dz -1)
    horizons[:,3]  = np.rint(pos_sal1[:,4]/dz -1)
    # creating horizon between fundomar and pos_sal1
    horizons[:,2] = np.rint((horizons[:,1] + horizons[:,3])/2)

    horizons[:,4]  = np.rint(pos_sal2[:,4]/dz -1)
    horizons[:,6]  = np.rint(pos_sal3[:,4]/dz -1)
    # creating horizon between pos_sal2 and pos_sal3
    horizons[:,5] = np.rint((horizons[:,4] + horizons[:,6])/2)

    horizons[:,9]  = np.rint(toposal[:,4]/dz -1)
    # creating horizon between pos_sal3 and toposal
    horizons[:,7] = np.rint((horizons[:,6] + horizons[:,9])/2)
    # creating horizon between pos_sal3.1 and toposal
    horizons[:,8] = np.rint((horizons[:,7] + horizons[:,9])/2)

    horizons[:,10]  = np.rint(basesal[:,4]/dz -1)
    
    # repeat base of salt horizon for layer inside reservoir salt's base + 170m
    horizons[:,11]  = np.rint(qc.smooth1D(basesal[:,4],50)/dz -1) + 17

    # repeat base of salt horizon for base of reservoir at 370 m bellow 
    horizons[:,12]  = np.rint(qc.smooth1D(basesal[:,4],100)/dz -1) + 37

    ## set velocity for each layer  - tomo1   
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
    horizon_vel[10] = 5450
    horizon_vel[11] = 4500
    horizon_vel[12] = 4800
else:
    "No valid options. Stop the program."
    exit()

pl.figure()
ax = pl.gca()
im = pl.imshow(velocity,cmap= "jet",extent=[trace[0],trace[-1],depth[-1],depth[0]],aspect="equal")
pl.xlabel("distance (m)")
pl.ylabel("depth (m)")

#adjust color bar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
ax.xaxis.tick_top()
ax.xaxis.set_label_position("top")

# plot horizons
for layer in range(1,N_layers):        
        ax.plot(horizons[:,0]*dx,horizons[:,layer]*dz)     

cb = pl.colorbar(im, cax=cax)
cb.ax.invert_yaxis()
pl.show(block=False) 

getmeanvagarosity = True

if (getmeanvagarosity):
    ## get mean vagarosity from input model by layer 
    layer = 1 
    mean_vel = 0.0
    Nsamples = 0
    print("getting mean velocity layer", layer)
    for ii in range(0,Nx):
        for jj in range(0,Nz):        
            if (jj <= horizons[ii,layer]):                    
                tmp = 1/velocity[jj,ii] # convert velocity to vagarosity
                mean_vel = mean_vel + tmp
                Nsamples += 1

    horizon_vel[0]  = 1/(mean_vel/Nsamples) # convert vagarosity to velocity
    print(horizon_vel[0])

    # 2nd layer - (N-1)th layer
    for layer in range(2,N_layers):
        mean_vel = 0.0
        Nsamples = 0
        print("getting mean velocity layer ", layer)
        for ii in range(0,Nx):
            for jj in range(0,Nz):   
                #layer =2
                if (jj > horizons[ii,layer-1] and jj <= horizons[ii,layer]):                    
                    tmp = 1/velocity[jj,ii] # convert velocity to vagarosity
                    mean_vel = mean_vel + tmp
                    Nsamples += 1

        horizon_vel[layer-1]  = 1/(mean_vel/Nsamples) # convert vagarosity to velocity
        print(horizon_vel[layer-1])
    mean_vel = 0.0
    Nsamples = 0
    # last layer
    layer = N_layers
    print("getting mean velocity layer ", layer)
    for ii in range(0,Nx):
        for jj in range(0,Nz):
            if (jj > horizons[ii,layer-1]):
                tmp = 1/velocity[jj,ii] # convert velocity to vagarosity
                mean_vel = mean_vel + tmp
                Nsamples += 1

    horizon_vel[layer-1]  = 1/(mean_vel/Nsamples) # convert vagarosity to velocity
    print(horizon_vel[layer-1])
## Fill velocity model
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

pl.figure()
ax = pl.gca()
im = pl.imshow(velocity_horizon,cmap= "jet",extent=[trace[0],trace[-1],depth[-1],depth[0]],aspect="equal")
pl.xlabel("distance (m)")
pl.ylabel("depth (m)")

#adjust color bar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
ax.xaxis.tick_top()
ax.xaxis.set_label_position("top")

# plot horizons
for layer in range(1,N_layers):        
        ax.plot(horizons[:,0]*dx,horizons[:,layer]*dz)     

cb = pl.colorbar(im, cax=cax)
cb.ax.invert_yaxis()
pl.show(block=False) 

# write velocity P
qc.savebinaryfile(Nz,Nx,velocity_horizon,outfile+".bin")