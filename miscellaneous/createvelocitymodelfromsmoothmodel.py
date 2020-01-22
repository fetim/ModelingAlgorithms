#!/opt/anaconda3/bin/python
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
    # first layer
    layer = 0
    print("creating layer ", layer)
    for ii in range(0,dim2):
        for jj in range(0,dim1):                       
            if(inputmatrix[jj,ii] <= layer_value[layer] ):
                outputmatrix[jj,ii] = layer_value[layer]                
    # second up to N-1 layer
    for layer in range(1,N_layers):
        print("creating layer ", layer)
        for ii in range(0,dim2):
            for jj in range(0,dim1):             
                    if ((inputmatrix[jj,ii] >= layer_value[layer-1] ) and (inputmatrix[jj,ii] <= layer_value[layer])):
                        outputmatrix[jj,ii] = layer_value[layer]  
    # last layer
    layer = N_layers - 1    
    print("creating layer ", layer+1)
    for ii in range(0,dim2):
        for jj in range(0,dim1): 
            if(inputmatrix[jj,ii] >= layer_value[layer] ):
                outputmatrix[jj,ii] = layer_value[layer]                        
                        
    return outputmatrix                      

if __name__ =="__main__":
    '''
    Ambiente para teste     
    '''

    import qualitycontrolfunctions as qc

    # parameters
    name = "modelsBuzios/xline_3825_buzios_1001nz_2401nx_dz10m_dx25m.bin"
    Nx = 2401
    Nz = 1001
    dx = 25
    dz = 10
    
    trace = np.arange(0,dx*Nx,dx)    
    depth = np.arange(0,dz*Nz,dz)    
    
    # import velocity model
    velocity = qc.readbinaryfile(Nz,Nx,name)

    # Generate layered model using a smooth model
    N_layers = 50
    profile = velocity[:,int(Nx/2)]
    
    # create layered profile
    profile_layer = antismooth(profile,N_layers,"rms")
    
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
    
    #### generating layered velocity model    
    velocity_layered = generate_intervalar_model(velocity,layers_rms)    
    
    pl.imshow(velocity_layered)
    pl.show()

    # write velocity P
    outfile = name + '_layer_cake_from_smooth.bin'
    qc.savebinaryfile(Nz,Nx,velocity_layered,outfile)