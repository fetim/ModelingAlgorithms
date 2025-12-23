import cupy as cp
"""
CUDA Kernels for wavefield modeling
"""
laplacian_8th_source = r'''
extern "C" __global__
void laplacian_kernel_8th_raw(const float* u, float* lap, float dz, float dx, int Nz, int Nx) {
    
    // Calculate Global Index
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int total_size = Nz * Nx;

    if (i >= total_size) return;

    // 2D Coordinates
    int iz = i / Nx;
    int ix = i % Nx;

    // Boundary Check
    if (iz < 4 || iz >= Nz - 4 || ix < 4 || ix >= Nx - 4) {
        lap[i] = 0.0f;
        return;
    }

    // Coefficients
    const float c0 = -205.0f / 72.0f;
    const float c1 =    8.0f / 5.0f;
    const float c2 =   -1.0f / 5.0f;
    const float c3 =    8.0f / 315.0f;
    const float c4 =   -1.0f / 560.0f;

    // Vertical Derivative (Z-Direction) 
    float d2z = c1 * (u[i - Nx]   + u[i + Nx])   +
                c2 * (u[i - 2*Nx] + u[i + 2*Nx]) +
                c3 * (u[i - 3*Nx] + u[i + 3*Nx]) +
                c4 * (u[i - 4*Nx] + u[i + 4*Nx]) +
                c0 * u[i];
    
    d2z /= (dz * dz); 

    // Horizontal Derivative (X-Direction)
    float d2x = c1 * (u[i - 1] + u[i + 1]) +
                c2 * (u[i - 2] + u[i + 2]) +
                c3 * (u[i - 3] + u[i + 3]) +
                c4 * (u[i - 4] + u[i + 4]) +
                c0 * u[i];

    d2x /= (dx * dx);

    // 7. Sum results
    lap[i] = d2z + d2x;
}
'''

# Compile kernel once
laplacian_kernel_8th_raw = cp.RawKernel(laplacian_8th_source, 'laplacian_kernel_8th_raw')

absorb_code = r'''
extern "C" __global__
void DampingWavefieldKernel(float* Uf , float* Uc, int nb, int Nz, int Nx) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int total_size = Nz * Nx;
    if (i >= total_size) return;

    // Get 2D coordinates from linear index 'i'
    int iz = i / Nx;
    int ix = i % Nx;

    const float sb = 3.0f * nb;
    float fb = 0.0f;

    if (iz < nb) {
        fb = (nb - iz) / (sqrtf(2.0f) * sb);
        Uf[i] = Uf[i] * expf(-fb * fb);
        Uc[i] = Uc[i] * expf(-fb * fb);
    }
    else if (iz >= Nz - nb + 1) {
        fb = (iz - Nz + nb) / (sqrtf(2.0f) * sb);
        Uf[i] = Uf[i] * expf(-fb * fb);
        Uc[i] = Uc[i] * expf(-fb * fb);
    }
    else if (ix < nb) {
        fb = (nb - ix) / (sqrtf(2.0f) * sb);
        Uf[i] = Uf[i] * expf(-fb * fb);
        Uc[i] = Uc[i] * expf(-fb * fb);
    }
    else if (ix >= Nx - nb + 1) {
        fb = (ix - Nx + nb) / (sqrtf(2.0f ) * sb);
        Uf[i] = Uf[i] * expf(-fb * fb);
        Uc[i] = Uc[i] * expf(-fb * fb);
    }
}
'''
DampingWavefieldKernel = cp.RawKernel(absorb_code, 'DampingWavefieldKernel')