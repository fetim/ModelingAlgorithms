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

def function_interp_grid_bilinear(n_p,n_n,grid_p,grid_n,GRID_NODES_P,GRID_NODES_N,F):
    #numero do p
    # profundidade

    D = np.zeros([n_p,n_n])
    DP = GRID_NODES_P[2]-GRID_NODES_P[1]
    DN = GRID_NODES_N[2]-GRID_NODES_N[1]
    NP = len(GRID_NODES_P)
    NN = len(GRID_NODES_N)
    for p in range(0,n_p): 
        for n in range(0,n_n):
            ind_p = int(np.floor(grid_p[p]/DP))
            ind_n = int(np.floor(grid_n[n]/DN))

            if (ind_p<0):
                 ind_p = 1
            if (ind_n<0):
                 ind_n = 1
            if (ind_p>NP-2):
                 ind_p = NP-2
            if (ind_n>NN-2):
                 ind_n = NN-2

            P00 = F[ind_p,ind_n]
            P01 = F[ind_p,ind_n+1]
            P10 = F[ind_p+1,ind_n]
            P11 = F[ind_p+1,ind_n+1]
            P0 = GRID_NODES_P[ind_p]
            P1 = GRID_NODES_P[ind_p+1]
            N0 = GRID_NODES_N[ind_n]
            N1 = GRID_NODES_N[ind_n+1]

            D[p,n] = bilinear_interpolation_m(P00,P01,P10,P11,P0,P1,N0,N1,grid_p[p],grid_n[n])
    
    return D

if __name__ =="__main__":  

    # parameters
    folder = "/mnt/Desktop/Dados_Buzios/" 
    name = "xline_3824_buzios_1001z_1402x_dz10m_dx50m"
    Nx = 1202
    Nz = 1001
    dx = 50
    dz = 10

    # # new parameters
    dx_interp = 25
    dz_interp = 10

    # import velocity model
    velocity = qc.readbinaryfile(Nz,Nx,folder+name+".bin")
    qc.plotmatrix(velocity,'jet')

    depth  = np.arange(0,Nz*dz,dz)
    east =  np.arange(0,Nx*dx,dx)

    Nx_interp = Nx*2 - 1 
    Nz_interp = Nz

    depth_interp  = np.arange(0,Nz_interp*dz_interp,dz_interp)
    east_interp =  np.arange(0,(Nx_interp)*dx_interp,dx_interp)    

    velocity_interp = function_interp_grid_bilinear(Nz_interp,Nx_interp,depth_interp,east_interp,depth,east,velocity)

    
    
    qc.plotmatrix(velocity_interp,'jet')            

    # write velocity P
    outfile =" xline_3824_buzios_1001z_2403x_dz10m_dx25m_interp"
    qc.savebinaryfile(Nz_interp,Nx_interp,velocity_interp,folder+outfile +".bin")
