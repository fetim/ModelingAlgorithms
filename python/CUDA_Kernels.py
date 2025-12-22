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
