from numpy import arange, zeros, pi, sqrt, exp, shape, copy, any, percentile, roll
from numba import jit,njit,prange
import cupy as cp

from CUDA_Kernels import LaplacianKernel8thOrder
from CUDA_Kernels import DampingWavefieldKernel
#%% Numba JIT-based functions
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
def laplacian_serial(Uc,dz,dx):
    Nz,Nx = Uc.shape
    lap = zeros((Nz,Nx))
    dx2 = 12.0 * dx * dx
    dz2 = 12.0 * dz * dz
    for i in range(2,Nx-3):
        for j in range(2,Nz-3):
            pxx = (-Uc[j,i+2] + 16*Uc[j,i+1] -30*Uc[j,i] + 16*Uc[j,i-1] - Uc[j,i-2] )/(dx2)
            pzz = (-Uc[j+2,i] + 16*Uc[j+1,i] -30*Uc[j,i] + 16*Uc[j-1,i] - Uc[j-2,i] )/(dz2)
            lap[j,i]  = (pxx + pzz)
    return lap


@staticmethod
@jit(nopython=True,parallel=True)
def laplacian_numba(Uc,dz,dx):
    """
    4th-order centered finite differences using explicit loops (numba JIT)
    """
    Nz,Nx = shape(Uc)
    lap = zeros((Nz,Nx))
    dx2 = 12.0 * dx * dx
    dz2 = 12.0 * dz * dz
    for i in prange(2,Nx-3): # parallel
        for j in prange(2,Nz-3): # parallel
            pxx = (-Uc[j,i+2] + 16*Uc[j,i+1] -30*Uc[j,i] + 16*Uc[j,i-1] - Uc[j,i-2] )/(dx2)
            pzz = (-Uc[j+2,i] + 16*Uc[j+1,i] -30*Uc[j,i] + 16*Uc[j-1,i] - Uc[j-2,i] )/(dz2)
            lap[j,i]  = (pxx + pzz)
    
    return lap

@staticmethod
def LaplacianCUDA(Uc_g,lap,dz,dx):
    Nz, Nx = Uc_g.shape
    
    # Grid setup
    total_pixels = Nz * Nx
    threads_per_block = 256
    blocks_per_grid = (total_pixels + threads_per_block - 1) // threads_per_block
    
    # Launch Kernel
    LaplacianKernel8thOrder(
        (blocks_per_grid,),       # Grid dimensions
        (threads_per_block,),     # Block dimensions
        (
            Uc_g,                # const float* u
            lap,                  # float* lap
            cp.float32(dz),       # float dz
            cp.float32(dx),       # float dx
            cp.int32(Nz),         # int Nz
            cp.int32(Nx)          # int Nx
        )
    )
    return lap

@staticmethod
def laplacianPseudoSpectralCUDA(Uc_g,lap,dz,dx):
    Uk = cp.fft.fft2(Uc_g)
    lapk = Uk * lap
    return cp.fft.ifft2(lapk,s=Uc_g.shape).real

@staticmethod
def acousticWaveEquationCUDA(Uf_g,Uc_g,lap,vp_g,dz,dx,dt,nb):
    Nz, Nx = Uf_g.shape
    
    # Acoustic Wavefield Update
    lap = LaplacianCUDA(Uc_g,lap,dz,dx)
    Uf_g = (vp_g * vp_g) * (dt * dt) * lap + 2.0 * Uc_g - Uf_g
    Uf_g,Uc_g = apply_dampingCUDA(Uf_g,Uc_g,nb)
    
    # Uf_g = (vp_g * vp_g) * (dt * dt) * laplacianPseudoSpectralCUDA(Uc_g,lap,dz,dx) + 2.0 * Uc_g - Uf_g
    
    
    return Uf_g,Uc_g

def apply_dampingCUDA(Uf_g,Uc_g,nb):
    Nz,Nx = cp.int32(cp.shape(Uf_g))
    #block and grid sizes
    total_size = Nz * Nx
    treads_per_block = 256
    # how many blocks are needed to cover all elements
    blocks_per_grid = (total_size + treads_per_block - 1) // treads_per_block
    #launch kernel
    DampingWavefieldKernel((blocks_per_grid,), (treads_per_block,), (Uf_g,Uc_g,nb,Nz,Nx))
    return Uf_g,Uc_g
    
