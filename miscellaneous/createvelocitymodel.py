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
    Ambiente para teste das funções criadas      
    '''
    import qualitycontrolfunctions as qc
    


    inputname = "velocitymodel_Marmousi_Nx560_Nz281_dh50_smth50.bin"

    Nx = 560
    Nz = 281
    dx = 50
    dz = 50
    velocity = qc.readbinaryfile(Nz,Nx,inputname)
    qc.plotmatrix(velocity,'jet')

    #############################################
    Nlayers = 10
    inputdata = velocity[:,200]
    profile = antismooth(inputdata,Nlayers,"rms")

    velocity_layer = np.zeros([Nz,Nx])
    for ii in range(0,Nx):
        velocity_layer[:,ii] = profile[:]     
    qc.plotmatrix(velocity_layer,'jet')

    pl.figure()
    depth = np.arange(0,dz*Nz,dz)    
    pl.plot(depth,inputdata,depth,profile)  

    #############################################
    N_elements = int(np.size(inputdata)/Nlayers)    
    resto = np.size(inputdata)%Nlayers

    layers_max  = np.zeros([Nlayers,1])
    layers_min  = np.zeros([Nlayers,1])
    layers_mean = np.zeros([Nlayers,1])
    layers_rms  = np.zeros([Nlayers,1])
    for layer in range(0,Nlayers):    
        layers_max[layer]  = np.max(inputdata[layer*N_elements:(layer+1)*N_elements])
        layers_min[layer]  = np.min(inputdata[layer*N_elements:(layer+1)*N_elements])
        layers_mean[layer] = np.mean(inputdata[layer*N_elements:(layer+1)*N_elements])
        layers_rms[layer]  = np.sqrt(np.mean(inputdata[layer*N_elements:(layer+1)*N_elements]**2))
     
    pl.figure()
    pl.plot(layers_rms)
    #############################################
    velocity_layer_v2 = np.zeros([Nz,Nx])
    for ii in range(0,Nx):
        for jj in range(0,Nz):
            layer = 0
            if(velocity[jj,ii] <= layers_rms[layer] ):
                velocity_layer_v2[jj,ii] = layers_rms[layer]
                # print("ok in      " ,layer,jj, ii, velocity_layer_v2[jj,ii],layers_rms[layer] )                

            layer = Nlayers - 1
            if(velocity[jj,ii] >= layers_rms[layer] ):
                velocity_layer_v2[jj,ii] = layers_rms[layer]

            else:   
                for layer in range(1,Nlayers):
                    if ((velocity[jj,ii] >= layers_rms[layer-1] ) and (velocity[jj,ii] <= layers_rms[layer])):
                        velocity_layer_v2[jj,ii] = layers_rms[layer]
                        # print("ok in      " ,layer,jj, ii, velocity_layer_v2[jj,ii],layers_rms[layer] )
                    # else:
                        # print("problem in ", layer,jj, ii, velocity[jj,ii],layers_rms[layer] )

    qc.plotmatrix(velocity_layer_v2,'jet')

    pl.show()  