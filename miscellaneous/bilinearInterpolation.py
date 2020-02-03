#!/opt/anaconda3/bin/python
import matplotlib.pyplot as pl
import numpy as np
import qualitycontrolfunctions as qc
from mpl_toolkits.axes_grid1 import make_axes_locatable

def bilinear_interpolation_m(P00,P01,P10,P11,P0,P1,N0,N1,p,n):
    A = (N1-N0)*(P1-P0)
    p = P00*(P1-p)*(N1-n) + P01*(P1-p)*(n-N0) + P10*(p-P0)*(N1-n) + P11*(p-P0)*(n-N0)
    p = p/A
    return p   

if __name__ =="__main__":  

    # parameters
    folder = "modelsBuzios/" 
    name = folder + "xline_buzios_1001z_1202x.bin"
    Nx = 1202
    Nz = 1001
    dx = 50
    dz = 10

    # import velocity model
    velocity = qc.readbinaryfile(Nz,Nx,name)
    qc.plotmatrix(velocity,'jet')

    # new parameters
    dx_new = 25
    dz_new = 10

    ratio_dx = dx/dx_new
    ratio_dz = dz/dz_new

    Nx_new = int(ratio_dx*Nx)
    Nz_new = int(ratio_dz*Nz)

    velocity_interp = np.zeros([Nz_new,Nx_new])

    # copy input to expanded model
    for ii in range(0,Nx):
        for jj in range(0,Nz):
            velocity_interp[int(ratio_dz*jj),int(ratio_dx*ii)] = velocity[jj,ii]

    ## Filling with bilinear interpolation
    for ii in range(0,Nx-1):
        # interp last line
        jj = Nz-1
        value = bilinear_interpolation_m(velocity[jj-1,ii-1],velocity[jj-1,ii],\
                velocity[jj,ii-1],velocity[jj,ii],0,1,0,1,0.5,0.5)
        velocity_interp[int(ratio_dz*jj),int(ratio_dx*ii)+1] = value                
        for jj in range(0,Nz-1):            
            value = bilinear_interpolation_m(velocity[jj,ii],velocity[jj,ii+1],\
                velocity[jj+1,ii],velocity[jj+1,ii+1],0,1,0,1,0.5,0.5)

            velocity_interp[int(ratio_dz*jj),int(ratio_dx*ii)+1] = value

    
    # interp last column
    ii=Nx-1
    for jj in range(0,Nz):        
        value = bilinear_interpolation_m(velocity[jj-1,ii-1],velocity[jj-1,ii],\
                velocity[jj,ii-1],velocity[jj,ii],0,1,0,1,0.5,0.5)
        velocity_interp[int(ratio_dz*jj),int(ratio_dx*ii)+1] = value                

    qc.plotmatrix(velocity_interp,'jet')            

# write velocity P
outfile = folder + "xline_buzios_1001z_2401_interp.bin"
qc.savebinaryfile(Nz_new,2401,velocity_interp[:,0:2401],outfile)
