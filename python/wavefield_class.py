from numpy import arange, zeros, pi, sqrt, exp, shape, copy, any, percentile 
from numba import jit,njit,prange
import matplotlib.pyplot as plt
import time

class wavefield:
    approximation = "acoustic"

    def __init__(self):
        self.parameterpath = "parameter.txt"
        self.outputpath    =  ""

        self = self.readparameters()

    def readparameters(self):
        self.dx   = 5.
        self.dz   = 5.
        self.dt   = 0.001
 
        self.L    = 10000
        self.D    = 3000
        self.T    = 2

        self.fcut = 60 # Max frequency

        self.x    = arange(0,self.L+self.dx,self.dx)
        self.z    = arange(0,self.D+self.dz,self.dz)
        self.t    = arange(0,self.T+self.dt,self.dt)
 
        self.Nx   = len(self.x)
        self.Nz   = len(self.z)
        self.Nt   = len(self.t)

        self.Nrec = self.Nx

        self.vp         = zeros([self.Nz,self.Nx])
        self.current    = zeros([self.Nz,self.Nx])
        self.past       = zeros([self.Nz,self.Nx])
        self.future     = zeros([self.Nz,self.Nx])
        self.seismogram = zeros([self.Nt,self.Nrec])

        self.source     = self.wavelet_ricker(self.t,self.fcut)

        self.snap       = False

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

        

    def load_source(self):
        self.source = self.wavelet_ricker(self.t,60)        

    
    def acquisition_geometry(self):
        self.sx = int(self.L/2/self.dx)
        self.sz = 2

        self.rx = arange(0,self.Nrec,dtype=int)
        self.rz = zeros([self.Nrec],dtype=int) + 10

    def calculate_wavefield(self):
        self.future = self.update_equation(self.future,
                                           self.current,
                                           self.past,
                                           self.vp,
                                           self.dz,
                                           self.dx,
                                           self.dt)
        return self                                           

    def update_wavefield(self):
        self.past    = copy(self.current)
        self.current = copy(self.future)
        return self

    def register_seismogram(self,k):
        self.seismogram[k,self.rx] = self.current[self.rz,self.rx]
        return self
        
    def forward_modeling(self,sz,sx):        
        if self.snap: fig, ax = plt.subplots()
        for k in range(0,self.Nt):
            self.current[sz,sx] = self.current[sz,sx] - (self.dt*self.dt)*(self.vp[sz,sx]*self.vp[sz,sx]) * self.source[k]
            self.calculate_wavefield()
            self.update_wavefield()
            self.register_seismogram(k)
            if (k%100 ==0):
                print("step = %i" %k)
                if self.snap:
                    ax.cla()
                    ax.imshow(self.current)
                    plt.pause(0.001)           
        return self



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

    
    @staticmethod
    @jit(nopython=True,parallel=True)
    def velocitymodel_3layers(vp,v1,v2,v3):
        Nz,Nx = shape(vp)
        vp[0:int(Nz/4),:]         = v1
        vp[int(Nz/4):int(Nz/2),:] = v2
        vp[int(Nz/2):,:]          = v3
        return vp
        
    @staticmethod
    @jit(nopython=True,parallel=True)
    def wavelet_ricker(time,freq):        
        td     = time - 2*sqrt(pi)/freq
        fcd    = freq/(sqrt(pi)*3) 
        source = (1 - 2*pi*(pi*fcd*td)*(pi*fcd*td))*exp(-pi*(pi*fcd*td)*(pi*fcd*td))
        return source

    @staticmethod
    @jit(nopython=True,parallel=True)
    def update_equation(Uf,Uc,Up,vp,dz,dx,dt):
        Nz,Nx = shape(Uf)
        for i in prange(2,Nx-3): # parallel
            for j in prange(2,Nz-3): # parallel
                pxx = (-Uc[j,i+2] + 16*Uc[j,i+1] -30*Uc[j,i] + 16*Uc[j,i-1] - Uc[j,i-2] )/(12*dx*dx)
                pzz = (-Uc[j+2,i] + 16*Uc[j+1,i] -30*Uc[j,i] + 16*Uc[j-1,i] - Uc[j-2,i] )/(12*dz*dz)
                Uf[j,i]  = (vp[j,i]*vp[j,i])*(dt*dt)*(pxx + pzz) + 2*Uc[j,i] - Up[j,i]
        return Uf


start_time = time.time()

# Start parameters and wavefields
u = wavefield()

# load velocity model
u.vp = u.velocitymodel_3layers(u.vp,1500,2000,3000)

u.check_dispersionstability()

# generate synthetic seismogram
u.forward_modeling(u.sz,u.sx)

u.plot_seismogram()

print("Normal end of execution")
print("--- %s seconds ---" % (time.time() - start_time))

    



