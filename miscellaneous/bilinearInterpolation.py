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

    dim1=2
    dim2=2
    matrix =  np.zeros([dim1,dim2])

    matrix[0,0] = 2
    matrix[0,1] = 2
    matrix[1,0] = 4
    matrix[1,1] = 4

    interpolation = bilinear_interpolation_m(matrix[0,0],matrix[0,1],matrix[1,0],matrix[1,1],0,1,0,1,0.5,0.5)
    print(interpolation)

     # parameters
    name = "modelsBuzios/xline_buzios_1001z_1202x.bin"
    Nx = 1202
    Nz = 1001
    dx = 50
    dz = 10
    
    trace = np.arange(0,dx*Nx,dx)    
    depth = np.arange(0,dz*Nz,dz)    
    
    # import velocity model
    velocity = qc.readbinaryfile(Nz,Nx,name)

    qc.plotmatrix(velocity,'jet')

    dx_new = 25
    dz_new = 10

    ratio_dx = dx/dx_new
    ratio_dz = dz/dz_new

    Nx_new = int(ratio_dx*Nx)
    Nz_new = int(ratio_dz*Nz)

    trace_new = np.arange(0,dx_new*Nx,dx_new)    
    depth_new = np.arange(0,dz_new*Nz,dz_new)


    velocity_interp = np.zeros([Nz_new,Nx_new])

    # copy input to expanded model
    for ii in range(0,Nx):
        for jj in range(0,Nz):
            velocity_interp[int(ratio_dz*jj),int(ratio_dx*ii)] = velocity[jj,ii]




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

    ii=Nx-1
    # last column
    for jj in range(0,Nz):
        # print(velocity[jj,ii-1],velocity[jj,ii],\
        #         velocity[jj+1,ii-1],velocity[jj+1,ii])
        value = bilinear_interpolation_m(velocity[jj-1,ii-1],velocity[jj-1,ii],\
                velocity[jj,ii-1],velocity[jj,ii],0,1,0,1,0.5,0.5)
        velocity_interp[int(ratio_dz*jj),int(ratio_dx*ii)+1] = value                

    qc.plotmatrix(velocity_interp,'jet')            
