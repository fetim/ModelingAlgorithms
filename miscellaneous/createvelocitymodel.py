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

def generate_intervalar_model(inputmatrix,layer_value):
    '''
    Generate velocity model using rms velocity reference
    '''
    dim1, dim2 = np.shape(inputmatrix)
    outputmatrix = np.zeros([dim1,dim2])
    N_layer = np.shape(layer_value)
    for ii in range(0,dim2):
        for jj in range(0,dim1):
            # first layer
            layer = 0
            if(inputmatrix[jj,ii] <= layer_value[layer] ):
                outputmatrix[jj,ii] = layer_value[layer]                
            # last layer
            layer = N_layers - 1
            if(inputmatrix[jj,ii] >= layer_value[layer] ):
                outputmatrix[jj,ii] = layer_value[layer]
            # second up to N-1 layer
            else:   
                for layer in range(1,N_layers):
                    if ((inputmatrix[jj,ii] >= layer_value[layer-1] ) and (inputmatrix[jj,ii] <= layer_value[layer])):
                        outputmatrix[jj,ii] = layer_value[layer]  

    return outputmatrix                      

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

    # # Check velocity model
    # qc.plotmatrix(velocity,'jet')


    N_layers = 10
    profile = velocity[:,200]
    
    # create layered profile
    profile_layer = antismooth(profile,N_layers,"rms")
    
    # compare profiles
    # pl.figure()
    # depth = np.arange(0,dz*Nz,dz)    
    # pl.plot(depth,profile,depth,profile_layer)  

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
    
    # # check velocity layers
    pl.figure()
    pl.plot(layers_min,'*')

    # # generating layered velocity model
    velocity_layer = generate_intervalar_model(velocity,layers_min)    

    # # check layered velocity model
    qc.plotmatrix(velocity_layer,'jet')

    N_horizon = 1
    N_layers = N_horizon + 1
    horizon = np.zeros([Nx,2,N_horizon])

    aux = np.zeros([Nx,2])  
    x = np.arange(0,Nx)
    k = 2*np.pi/(Nx)
    phi = 5*np.pi/2
    z0 = 5*(np.sin(k * x - phi) + 1) + 50
    
    aux[:,0] = x
    aux[:,1] = z0
    horizon[:,:,0] = np.rint(aux)   
    
    # z1 = 5*(np.sin(k * x - 3/4*phi) + 1) + 90
    # aux[:,1] = z1
    # horizon[:,:,1] = np.rint(aux)

    # z2 = 5*(np.sin(k * x - 3*phi) + 1) + 150
    # aux[:,1] = z2
    # horizon[:,:,2] = np.rint(aux)  
    
    # Edit first layer of the model
    if(N_horizon == 1):
        for ii in range(0,Nx):
            for jj in range(0,Nz):
                # first layer
                layer = 0 
                if (jj <= horizon[ii,1,0]):
                    velocity_layer[jj,ii] = layers_min[layer]

    # Create velocity model using layers
    elif (N_horizon >= 2):
        for ii in range(0,Nx):
            for jj in range(0,Nz):
                # first layer
                layer = 0 
                if (jj <= horizon[ii,1,0]):
                    velocity_layer[jj,ii] = layers_min[layer]
                # Last layer
                layer = N_layers - 1 
                if (jj > horizon[ii,1,N_horizon-1]):
                    velocity_layer[jj,ii] = layers_min[layer]
                # second up to N-1 layer    
                else:
                    for layer in range(1,N_layers-1):                    
                        if ((jj > horizon[ii,1,layer-1]) and jj <= horizon[ii,1,layer]):                        
                            velocity_layer[jj,ii] = layers_min[layer]


    # check layered velocity model
    qc.plotmatrix(velocity_layer,'jet')    
    for layer in range(0,N_layers-1):
        pl.plot(horizon[:,0,layer],horizon[:,1,layer])
    
    pl.show()  