from numpy import arange, zeros, pi, sqrt, exp, shape, copy, any, percentile, roll
from numba import jit,njit,prange
import cupy as cp

from CUDA_Kernels import laplacian_kernel
from CUDA_Kernels import laplacian_kernel_8th   
from CUDA_Kernels import laplacian_kernel_8th_raw
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
def laplacian_roll(Uc,dz,dx):
    """
    4th-order centered finite differences using roll (vectorized)
    """
    Nz,Nx = shape(Uc)
    dx2 = 12.0 * dx * dx
    dz2 = 12.0 * dz * dz
    # rolls wrap; we'll zero the outer 2 rows/cols after computation to avoid wrap artifacts
    pxx = (-roll(Uc, -2, axis=1) + 16.0 * roll(Uc, -1, axis=1)
            - 30.0 * Uc + 16.0 * roll(Uc, 1, axis=1) - roll(Uc, 2, axis=1)) / dx2
    
    pzz = (-roll(Uc, -2, axis=0) + 16.0 * roll(Uc, -1, axis=0)
            - 30.0 * Uc + 16.0 * roll(Uc, 1, axis=0) - roll(Uc, 2, axis=0)) / dz2
    
    lap = pxx + pzz
    # zero the borders (2 layers) to avoid wrap contributions
    lap[:2, :] = 0
    lap[-2:, :] = 0
    lap[:, :2] = 0
    lap[:, -2:] = 0
    return (lap)

@staticmethod
def laplacian_cupy(Uc,dz,dx):
    """
    4th-order centered finite differences using CuPy
    """
    Nz,Nx = shape(Uc)
    Uc_g = cp.asarray(Uc)
    lap = cp.zeros_like(Uc_g)
    denom_x = 12.0 * dx * dx
    denom_z = 12.0 * dz * dz    
    
    pxx = (-cp.roll(Uc_g, -2, axis=1) + 16.0 * cp.roll(Uc_g, -1, axis=1)
            - 30.0 * Uc_g + 16.0 * cp.roll(Uc_g, 1, axis=1) - cp.roll(Uc_g, 2, axis=1)) / denom_x
    pzz = (-cp.roll(Uc_g, -2, axis=0) + 16.0 * cp.roll(Uc_g, -1, axis=0)
            - 30.0 * Uc_g + 16.0 * cp.roll(Uc_g, 1, axis=0) - cp.roll(Uc_g, 2, axis=0)) / denom_z

    lap = pxx + pzz
    # Uf_g[...] = (vp_g * vp_g) * (dt * dt) * (pxx + pzz) + 2.0 * Uc_g - Up_g
    
    # zero the borders (2 layers) to avoid wrap contributions
    lap[:2, :] = 0
    lap[-2:, :] = 0
    lap[:, :2] = 0
    lap[:, -2:] = 0
    return cp.asnumpy(lap)

#%% CUDA kernel-based update function
def acousticWaveEquationCUDA(Uf_g,Uc_g,vp_g,dz,dx,dt):
    Nz,Nx = cp.int32(cp.shape(Uf_g))
    lap = cp.zeros_like(Uf_g)
    # laplacian_kernel(Uc_g, dz, dx, Nz, Nx, lap)
    laplacian_kernel_8th(Uc_g, dz, dx, Nz, Nx, lap)
    Uf_g = (vp_g * vp_g) * (dt * dt) * lap + 2.0 * Uc_g - Uf_g
    return Uf_g

def acousticWaveEquationCUDA_raw(Uf_g,Uc_g,vp_g,dz,dx,dt):
    Nz, Nx = Uf_g.shape
    
    # Pre-allocate output
    lap = cp.zeros_like(Uc_g)
    
    # Grid setup
    total_pixels = Nz * Nx
    threads_per_block = 256
    blocks_per_grid = (total_pixels + threads_per_block - 1) // threads_per_block
    
    # Launch Kernel
    laplacian_kernel_8th_raw(
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
    return (vp_g * vp_g) * (dt * dt) * lap + 2.0 * Uc_g - Uf_g


def update_wavefieldCUDA(Uc_g,Uf_g):
    Uf_g,Uc_g = cp.copy(Uc_g), cp.copy(Uf_g)
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
    
