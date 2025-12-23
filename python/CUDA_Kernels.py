import cupy as cp
"""
CUDA Kernels for wavefield modeling
"""

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

laplacian_kernel_8th = cp.ElementwiseKernel(
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

    const float c0 = -205.0f / 72.0f;
    const float c1 = 8.0f / 5.0f;
    const float c2 = -1.0f / 5.0f;
    const float c3 = 8.0f / 315.0f;
    const float c4 = -1.0f / 560.0f;

    T dfdz = c1 * (u[i-1] + u[i+1]) +
             c2 * (u[i-2] + u[i+2]) +
             c3 * (u[i-3] + u[i+3]) +
             c4 * (u[i-4] + u[i+4]) +
            c0 * u[i];
    dfdz /= (dz * dz);

    T dfdx = c1 * (u[i-Nx] + u[i+Nx]) +
             c2 * (u[i-2*Nx] + u[i+2*Nx]) +
             c3 * (u[i-3*Nx] + u[i+3*Nx]) +
             c4 * (u[i-4*Nx] + u[i+4*Nx]) +
            c0 * u[i];
    dfdx /= (dx * dx);

    lap = dfdz + dfdx;
    ''',
    'laplacian_kernel'
    )


# ABKernel = cp.ElementwiseKernel(
#     # 1. Inputs: We use 'raw' to modify arrays in-place
#     'raw T Uf, raw T Uc, int32 nb, int32 Nz, int32 Nx',
    
#     # 2. Outputs: Empty string (we are not returning a new array)
#     'T Ufab, T Ucab',
    
#     # 3. C++ Operation
#     '''
#     // Get 2D coordinates from linear index 'i'
#     int iz = i / Nx;
#     int ix = i % Nx;

#     // Safety check (should not happen if size is set correctly, but good practice)
#     if (iz < 0 || iz >= Nz || ix < 0 || ix >= Nx) {
#         return;
#     }

#     const float sb = 3.0f * nb;
#     float fb = 0.0f;

#     if (iz < nb) {
#         fb = (nb - iz) / (sqrtf(2.0f) * sb);
#         Ufab[i] = Uf[i] * expf(-fb * fb);
#         Ucab[i] = Uc[i] * expf(-fb * fb);
#     } 
#     ''',
#     'ABKernel'
# )

absorb_code = r'''
extern "C" __global__
void ABKernel(float* Uf , float* Uc, int nb, int Nz, int Nx) {
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
ABKernel = cp.RawKernel(absorb_code, 'ABKernel')