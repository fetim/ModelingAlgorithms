from numpy import arange, zeros, pi, sqrt, exp, shape, copy, any, percentile, roll
from numba import jit,njit,prange
import matplotlib.pyplot as plt
import time
import cupy as cp
import numpy as np

# from CUDA_Kernels import laplacian_kernel
from solvers import velocitymodel_3layers
from solvers import wavelet_ricker
from solvers import laplacian_serial
from solvers import laplacian_numba
from solvers import laplacian_roll
from solvers import laplacian_cupy

from solvers import acousticWaveEquationCUDA_raw
from solvers import apply_dampingCUDA
class wavefield:

    def __init__(self):
        self.parameterpath = "parameter.txt"
        self.outputpath    =  ""

        self = self.readparameters()

    def readparameters(self):
        self.dx   = 5.
        self.dz   = 5.
        self.dt   = 0.0005
 
        self.L    = 5000
        self.D    = 2000
        self.T    = 1

        self.fcut = 60 # Max frequency

        self.x    = arange(0,self.L+self.dx,self.dx)
        self.z    = arange(0,self.D+self.dz,self.dz)
        self.t    = arange(0,self.T+self.dt,self.dt)
 
        self.Nx   = len(self.x)
        self.Nz   = len(self.z)
        self.Nt   = len(self.t)

        self.nb   = 50  # absorbing boundary thickness in grid points
        self.Nrec = self.Nx

        self.vp         = zeros([self.Nz,self.Nx])
        
        self.source     = wavelet_ricker(self.t,self.fcut)

        self.snap       = True

        self.acquisition_geometry()      
                
        print("wavefield started")
        print("Nx = %i; Nz=%i; Nt =%i" %(self.Nx,self.Nz,self.Nt))
        return self

    
    def check_dispersionstability(self):
        '''' 2nd order in time and 4th order in space'''
        pointsWavelength = 5
        a_t              = 1 + 2 + 1 
        a_s              = 1/12 + 4/3 + 5/2 + 4/3 + 1/12

        lamb_critical    = self.vp.min()/self.fcut
        self.dx_critical = lamb_critical/pointsWavelength        
        self.dt_critical = self.dx_critical * sqrt(a_t/a_s)/self.vp.max()

        print("Courant number = %f" %(self.dt*self.vp.max()/self.dx))
        print("sqrt(a_t/a_s)  = %f" %sqrt(a_t/a_s))


        if (self.dt<=self.dt_critical and self.dx<=self.dx_critical):
            print("Stability and non-dispersion conditions ...... OK!") 
        else:
            print("Problem!")
            print("dt_critical = %f dt = %f" %(self.dt_critical,self.dt))
            print("dx_critical = %f dx = %f" %(self.dx_critical,self.dx))

    def acquisition_geometry(self):
        self.sx = int(self.L/2/self.dx)
        self.sz = 5

        self.rx = arange(0,self.Nrec,dtype=int)
        self.rz = zeros([self.Nrec],dtype=int) + 10
        self.Nrec = len(self.rx)
        self.channel = np.arange(0,self.Nrec,1)

#%% Memory allocation
    def allocate_wavefields(self):
        self.current    = zeros([self.Nz,self.Nx])
        self.future     = zeros([self.Nz,self.Nx])  
        self.seismogram = zeros([self.Nt,self.Nrec])

#%% CPU-based forward modeling
    def acousticWaveEquation(self):
        """
        Solve acoustic wave equation using 4th-order centered finite differences
        """
        # lap = laplacian_serial(self.current,self.dz,self.dx)
        # lap = laplacian_numba(self.current,self.dz,self.dx)
        # lap = laplacian_roll(self.current,self.dz,self.dx)
        lap = laplacian_cupy(self.current,self.dz,self.dx)
        
        self.future = (self.vp * self.vp) * (self.dt * self.dt) * lap + 2*self.current - self.past
                                                   
    def update_wavefield(self):
        self.past    = copy(self.current)
        self.current = copy(self.future)
    
    def register_seismogram(self,k):
        self.seismogram[k,self.rx] = self.current[self.rz,self.rx]
        
    def forward_modeling(self,sz,sx):       
        print("info: CPU Forward Modeling")
        if self.snap: fig, ax = plt.subplots()
        
        self.allocate_wavefields()        
        for k in range(0,self.Nt):
            self.current[sz,sx] = self.current[sz,sx] - (self.dt*self.dt)*(self.vp[sz,sx]*self.vp[sz,sx]) * self.source[k]
            self.acousticWaveEquation()
            self.update_wavefield()
            self.register_seismogram(k)
            if (k%100 ==0):
                print("step = %i" %k)
                if self.snap:
                    ax.cla()
                    ax.imshow(self.current)
                    ax.plot(sx,sz,'r+')
                    plt.pause(0.001)           
        return self.seismogram
#%% GPU-based forward modeling        
    def register_seismogramCUDA(self,k,Uc_g):
        self.seismogram[k,self.channel] = cp.asnumpy(Uc_g[self.rz,self.rx])

    def expandModel(self):

        nb = self.nb
        Nx,Nz = self.Nx,self.Nz

        self.Nx += 2*self.nb
        self.Nz += 2*self.nb
        self.sx += self.nb
        self.sz += self.nb
        self.rx += self.nb
        self.rz += self.nb
        
        vp_exp = np.zeros((Nz+2*nb,Nx+2*nb),dtype=np.float32)
        vp_exp[nb:-nb,nb:-nb] = self.vp
        vp_exp[0:nb,nb:-nb] = self.vp[0,:]
        vp_exp[-nb:,nb:-nb] = self.vp[-1,:]
        vp_exp[:,0:nb] = vp_exp[:,nb:nb+1]
        vp_exp[:,-nb:] = vp_exp[:,-nb-1:-nb]
        self.vp = vp_exp

    def forward_modelingGPU(self):
        print("info: GPU Forward Modeling")
        if self.snap: fig, ax = plt.subplots()

        self.expandModel()        
        self.allocate_wavefields()

        Uf_g = cp.zeros_like(cp.float32(self.future))
        Uc_g = cp.zeros_like(cp.float32(self.current))
        vp_g = cp.asarray(cp.float32(self.vp))
        source_g = cp.asarray(cp.float32(self.source))
        
        sz,sx = self.sz,self.sx
        for k in range(0,self.Nt):
            Uc_g[sz,sx] = Uc_g[sz,sx] - (self.dt*self.dt)*(vp_g[sz,sx]*vp_g[sz,sx]) * source_g[k]
            Uf_g = acousticWaveEquationCUDA_raw(Uf_g,Uc_g,vp_g,self.dz,self.dx,self.dt)
            Uf_g,Uc_g = Uc_g,Uf_g # swap pointers
            Uf_g,Uc_g = apply_dampingCUDA(Uf_g,Uc_g,self.nb)

            self.register_seismogramCUDA(k,Uc_g)

            if (k%100 ==0):
                print("step = %i" %k)
                if self.snap:
                    ax.cla()
                    ax.imshow(Uc_g[self.nb:-self.nb,self.nb:-self.nb].get())
                    plt.pause(0.001)           
        return self.seismogram

#%% Plotting functions
    def plot_velocitymodel(self):
        if any(self.vp):       
            plt.figure(figsize=(10,10))
            plt.imshow(self.vp,aspect="auto")
            plt.title('acoustic velocity')
            plt.show()
        else:
            print("Nothing to plot")

    def plot_sourcewavelet(self):
        if any(self.source):
            plt.figure()
            plt.plot(self.t,self.source)
            plt.show()
        else:
            print("Nothing to plot")

    def plot_seismogram(self):
        plt.figure(figsize=(10,10))
        perc = percentile(self.seismogram,99)
        plt.imshow(self.seismogram,aspect='auto',cmap='gray',vmin=-perc,vmax=perc)
        plt.show(block=False)

#%% Main execution
if __name__ == "__main__":
    start_time = time.time()

    # Start parameters and wavefields
    u = wavefield()

    # load velocity model
    u.vp = velocitymodel_3layers(u.vp,1500,2000,3000)
    u.check_dispersionstability()

    # generate synthetic seismogram
    # seismogram = u.forward_modeling(u.sz,u.sx)
    seismogram = u.forward_modelingGPU()

    u.plot_seismogram()

    print("Normal end of execution")
    print("--- %s seconds ---" % (time.time() - start_time))

    



