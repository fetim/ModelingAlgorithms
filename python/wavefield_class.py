from numpy import arange, zeros, pi, sqrt, exp, shape, copy, any, percentile, roll
from numba import jit,njit,prange
import matplotlib.pyplot as plt
import time
import cupy as cp


laplacian_kernel = cp.ElementwiseKernel(
    'raw T u, float32 dz, float32 dx, int32 Nz, int32 Nx',
    'T lap',
    '''
    // Calculate 2D coordinates from flat index
    int iz = i / Nx;
    int ix = i % Nx;

    // check boundaries 
    if (iz < 2 || iz >= Nz - 2 || ix < 2 || ix >= Nx - 2) {
        lap = 0.0f;
        return;
    }

    // 4th-order centered finite differences coefficients
    // [-1, 16, -30, 16, -1] / (12 * d*d)
    T dfdz = 16.0f * (u[i-1] + u[i+1]) -
              1.0f * (u[i-2] + u[i+2]) -
             30.0f * u[i];
    dfdz /= (12.0f * dz * dz);

    T dfdx = 16.0f * (u[i-Nx] + u[i+Nx]) -
              1.0f * (u[i-2*Nx] + u[i+2*Nx]) -
                30.0f * u[i];
    dfdx /= (12.0f * dx * dx);
    lap = dfdz + dfdx;
    ''',
    'laplacian_kernel'
    )

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

        

    def load_source(self):
        self.source = self.wavelet_ricker(self.t,60)        

    
    def acquisition_geometry(self):
        self.sx = int(self.L/2/self.dx)
        self.sz = 2

        self.rx = arange(0,self.Nrec,dtype=int)
        self.rz = zeros([self.Nrec],dtype=int) + 10

    def calculate_wavefield(self):
        self.future = self.update_equation_kernel(self.future,
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
    
    def calculate_wavefieldGPU(self,Uf_g,Uc_g,Up_g,vp_g,dz,dx,dt):
        Nz,Nx = shape(Uf_g)
        Uf_g = self.update_equation_kernelGPU(Uf_g,
                                             Uc_g,
                                             Up_g,
                                             vp_g,
                                             dz,
                                             dx,
                                             dt)
        return Uf_g
    
    def update_wavefieldGPU(self,Up_g,Uc_g,Uf_g):
        Up_g = cp.copy(Uc_g)
        Uc_g = cp.copy(Uf_g)
        return Up_g,Uc_g
    
    def register_seismogramGPU(self,k,Uc_g):
        self.seismogram[k,self.rx] = cp.asnumpy(Uc_g[self.rz,self.rx])


    def forward_modelingGPU(self,sz,sx):        
        if self.snap: fig, ax = plt.subplots()

        Uf_g = cp.zeros_like(cp.float32(self.future))
        Uc_g = cp.zeros_like(cp.float32(self.current))
        Up_g = cp.zeros_like(cp.float32(self.past))
        vp_g = cp.asarray(cp.float32(self.vp))
        source_g = cp.asarray(cp.float32(self.source))

        for k in range(0,self.Nt):
            Uc_g[sz,sx] = Uc_g[sz,sx] - (self.dt*self.dt)*(vp_g[sz,sx]*vp_g[sz,sx]) * source_g[k]
            Uf_g = self.calculate_wavefieldGPU(Uf_g,Uc_g,Up_g,vp_g,self.dz,self.dx,self.dt)

            
            Up_g,Uc_g = self.update_wavefieldGPU(Up_g,Uc_g,Uf_g)

            # self.calculate_wavefield)
            self.register_seismogramGPU(k,Uc_g)

            if (k%100 ==0):
                print("step = %i" %k)
                if self.snap:
                    ax.cla()
                    ax.imshow(Uc_g.get())
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
    
    @staticmethod
    def update_equation_roll(Uf,Uc,Up,vp,dz,dx,dt):
        # # Placeholder for potential alternative implementation
        Nz,Nx = shape(Uf)
        denom_x = 12.0 * dx * dx
        denom_z = 12.0 * dz * dz
        # rolls wrap; we'll zero the outer 2 rows/cols after computation to avoid wrap artifacts
        pxx = (-roll(Uc, -2, axis=1) + 16.0 * roll(Uc, -1, axis=1)
               - 30.0 * Uc + 16.0 * roll(Uc, 1, axis=1) - roll(Uc, 2, axis=1)) / denom_x
        
        pzz = (-roll(Uc, -2, axis=0) + 16.0 * roll(Uc, -1, axis=0)
               - 30.0 * Uc + 16.0 * roll(Uc, 1, axis=0) - roll(Uc, 2, axis=0)) / denom_z
        
        Uf[...] = (vp * vp) * (dt * dt) * (pxx + pzz) + 2.0 * Uc - Up
        return Uf
    
    @staticmethod
    def update_equation_cupy(Uf,Uc,Up,vp,dz,dx,dt):
        Nz,Nx = shape(Uf)
        Uf_g = cp.asarray(Uf)
        Uc_g = cp.asarray(Uc)
        Up_g = cp.asarray(Up)
        vp_g = cp.asarray(vp)
        denom_x = 12.0 * dx * dx
        denom_z = 12.0 * dz * dz    

        pxx = (-cp.roll(Uc_g, -2, axis=1) + 16.0 * cp.roll(Uc_g, -1, axis=1)
               - 30.0 * Uc_g + 16.0 * cp.roll(Uc_g, 1, axis=1) - cp.roll(Uc_g, 2, axis=1)) / denom_x
        pzz = (-cp.roll(Uc_g, -2, axis=0) + 16.0 * cp.roll(Uc_g, -1, axis=0)
               - 30.0 * Uc_g + 16.0 * cp.roll(Uc_g, 1, axis=0) - cp.roll(Uc_g, 2, axis=0)) / denom_z

        Uf_g[...] = (vp_g * vp_g) * (dt * dt) * (pxx + pzz) + 2.0 * Uc_g - Up_g
        
        # zero the borders (2 layers) to avoid wrap contributions
        Uf_g[:2, :] = 0
        Uf_g[-2:, :] = 0
        Uf_g[:, :2] = 0
        Uf_g[:, -2:] = 0
        return cp.asnumpy(Uf_g)
    
    

    @staticmethod
    def update_equation_kernel(Uf,Uc,Up,vp,dz,dx,dt):
        
        Nz,Nx = cp.int32(shape(Uf))
        Uf_g = cp.empty_like(cp.float32(Uf))
        Uc_g = cp.asarray(cp.float32(Uc))
        Up_g = cp.asarray(cp.float32(Up))
        vp_g = cp.asarray(cp.float32(vp))

        laplacian_kernel(Uc_g, dz, dx, Nz, Nx, Uf_g)
        Uf_g = (vp_g * vp_g) * (dt * dt) * Uf_g + 2.0 * Uc_g - Up_g
        
        return cp.asnumpy(Uf_g)
        
    @staticmethod
    def update_equation_kernelGPU(Uf_g,Uc_g,Up_g,vp_g,dz,dx,dt):
        
        Nz,Nx = cp.int32(cp.shape(Uf_g))
        laplacian_kernel(Uc_g, dz, dx, Nz, Nx, Uf_g)
        Uf_g = (vp_g * vp_g) * (dt * dt) * Uf_g + 2.0 * Uc_g - Up_g
        return Uf_g
    
if __name__ == "__main__":
    start_time = time.time()

    # Start parameters and wavefields
    u = wavefield()

    # load velocity model
    u.vp = u.velocitymodel_3layers(u.vp,1500,2000,3000)

    u.check_dispersionstability()

    # generate synthetic seismogram
    # u.forward_modeling(u.sz,u.sx)
    # 
    u.forward_modelingGPU(u.sz,u.sx)

    u.plot_seismogram()

    print("Normal end of execution")
    print("--- %s seconds ---" % (time.time() - start_time))

    



