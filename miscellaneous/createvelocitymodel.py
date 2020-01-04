#!/usr/bin/python
import matplotlib.pyplot as pl
import numpy as np

def antismooth(inputdata,N_layer,method):
    "Convert a 1D smooth function in a layered function"
    N_elements = int(np.size(inputdata)/N_layer)    
    resto = np.size(inputdata)%N_layer
    output = np.zeros(np.size(inputdata))
    
    if (method == "mean"):
        for layer in range(0,N_layer):    
            value = np.mean(inputdata[layer*N_elements:(layer+1)*N_elements])
            output[layer*N_elements:(layer+1)*N_elements] = value
        
        if (resto != 0):
            output[(layer+1)*N_elements:] = value

    elif(method == "rms"):
        for layer in range(0,N_layer):                
            value = np.sqrt(np.mean(inputdata[layer*N_elements:(layer+1)*N_elements]**2))
            output[layer*N_elements:(layer+1)*N_elements] = value
        
        if (resto != 0):
            output[(layer+1)*N_elements:] = value

    else:
        print("Using mean by default")
        for layer in range(0,N_layer):    
            value = np.mean(inputdata[layer*N_elements:(layer+1)*N_elements])
            output[layer*N_elements:(layer+1)*N_elements] = value
        
        if (resto != 0):
            output[(layer+1)*N_elements:] = value

    return output

if __name__ =="__main__":
    '''
    Ambiente para teste     
    '''
    import qualitycontrolfunctions as qc

    # parameters
    name = "velocitymodel_Marmousi_Nx560_Nz281_dh50_smth50.bin"
    Nx = 560
    Nz = 281
    dx = 50
    dz = 50
    velocity = qc.readbinaryfile(Nz,Nx,name)

    # Check velocity model
    qc.plotmatrix(velocity,'jet')


    N_layers = 7
    profile = velocity[:,200]
    
    # create layered profile
    profile_layer = antismooth(profile,N_layers,"rms")
    
    # compare profiles
    pl.figure()
    depth = np.arange(0,dz*Nz,dz)    
    pl.plot(depth,profile,depth,profile_layer)  

    # number of elements per layer 
    N_elements = int(np.size(profile)/N_layers)    
    resto = np.size(profile)%N_layers

    # generating velocity reference array
    layers_max  = np.zeros([N_layers,1])
    layers_min  = np.zeros([N_layers,1])
    layers_mean = np.zeros([N_layers,1])
    layers_rms  = np.zeros([N_layers,1])
    for layer in range(0,N_layers):    
        layers_max[layer]  = np.max(profile[layer*N_elements:(layer+1)*N_elements])
        layers_min[layer]  = np.min(profile[layer*N_elements:(layer+1)*N_elements])
        layers_mean[layer] = np.mean(profile[layer*N_elements:(layer+1)*N_elements])
        layers_rms[layer]  = np.sqrt(np.mean(profile[layer*N_elements:(layer+1)*N_elements]**2))
    
    # check velocity layers
    pl.figure()
    pl.plot(layers_rms,'*')

    
    # Generate velocity model using rms velocity reference
    velocity_layer = np.zeros([Nz,Nx])
    for ii in range(0,Nx):
        for jj in range(0,Nz):
            # first layer
            layer = 0
            if(velocity[jj,ii] <= layers_rms[layer] ):
                velocity_layer[jj,ii] = layers_rms[layer]                
            # last layer
            layer = N_layers - 1
            if(velocity[jj,ii] >= layers_rms[layer] ):
                velocity_layer[jj,ii] = layers_rms[layer]
            # second up to N-1 layer
            else:   
                for layer in range(1,N_layers):
                    if ((velocity[jj,ii] >= layers_rms[layer-1] ) and (velocity[jj,ii] <= layers_rms[layer])):
                        velocity_layer[jj,ii] = layers_rms[layer]                        

    # check velocity model layered
    qc.plotmatrix(velocity_layer,'jet')

    pl.show()  