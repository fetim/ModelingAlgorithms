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
    print("antes generetate")
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
    name = "xline_3825_buzios_1001nz_2401nx_dz10m_dx25m.bin"
    Nx = 2401
    Nz = 1001
    dx = 25
    dz = 10
    
    trace = np.arange(0,dx*Nx,dx)    
    depth = np.arange(0,dz*Nz,dz)    
    
    velocity = qc.readbinaryfile(Nz,Nx,name)

    # Check velocity model
    qc.plotmatrix(velocity,'jet')


    N_layers = 12
    profile = velocity[:,int(Nx/2)]
    
    # create layered profile
    profile_layer = antismooth(profile,N_layers,"rms")
    
    # # compare profiles
    # pl.figure()
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
    
    # # # check velocity layers
    # pl.figure()
    # pl.plot(layers_min,'*')
    # pl.show()

    # generating layered velocity model    
    velocity_layer = generate_intervalar_model(velocity,layers_rms)    
    
    # # check layered velocity model
    # qc.plotmatrix(velocity_layer,'jet')
    
    #load horizons
    fundomar = np.loadtxt("2D_horizons/Fundo_do_Mar_XL3825_2D_davidff.txt")
    pos_sal1 = np.loadtxt("2D_horizons/Pos-Sal_1_XL3825_2D_davidff.txt")
    pos_sal2 = np.loadtxt("2D_horizons/Pos-Sal_2_XL3825_2D_davidff.txt")
    pos_sal3 = np.loadtxt("2D_horizons/Pos-Sal_3_XL3825_2D_davidff.txt")
    toposal  = np.loadtxt("2D_horizons/Topo_Sal_XL3825_2D_davidff.txt")
    basesal  = np.loadtxt("2D_horizons/Base_Sal_XL3825_2D_davidff.txt")

    axis_x =  np.arange(0,Nx)
    fundomar_index = np.zeros([Nx,2])
    pos_sal1_index = np.zeros([Nx,2])
    pos_sal2_index = np.zeros([Nx,2])
    pos_sal3_index = np.zeros([Nx,2])
    toposal_index  = np.zeros([Nx,2])
    basesal_index  = np.zeros([Nx,2])

    fundomar_index[:,0] = axis_x
    pos_sal1_index[:,0] = axis_x
    pos_sal2_index[:,0] = axis_x
    pos_sal3_index[:,0] = axis_x
    toposal_index[:,0] = axis_x
    basesal_index[:,0] = axis_x

    # Check velocity model
    qc.plotmatrix(velocity,'jet')
    N_horizon = 6
    N_layers = N_horizon + 1
    horizons = np.zeros([Nx,N_layers])


    # Create horizon table
    horizons[:,0]  = np.arange(0,Nx)
    horizons[:,1]  = np.rint(fundomar[:-2,4]/dz -1)
    horizons[:,2]  = np.rint(pos_sal1[:-1,4]/dz -1)
    horizons[:,3]  = np.rint(pos_sal2[:,4]/dz -1)
    horizons[:,4]  = np.rint(pos_sal3[:-1,4]/dz -1)
    horizons[:,5]  = np.rint(toposal[:,4]/dz -1)
    horizons[:,6]  = np.rint(basesal[:,4]/dz -1)

    pl.plot(horizons[:,0],-horizons[:,1],'r')
    pl.plot(horizons[:,0],-horizons[:,2],'r')
    pl.plot(horizons[:,0],-horizons[:,3],'r')
    pl.plot(horizons[:,0],-horizons[:,4],'r')
    pl.plot(horizons[:,0],-horizons[:,5],'r')
    pl.plot(horizons[:,0],-horizons[:,6],'r')

    pl.show()  
        
    layers_value = np.zeros([N_layers,1])
    layers_value[0] = 1500
    layers_value[1] = 2500
    layers_value[2] = 3500
    layers_value[3] = 4500
    layers_value[4] = 5500
    layers_value[5] = 6500
    layers_value[6] = 7500

    # Edit first layer of the model
    for ii in range(0,Nx):
        for jj in range(0,Nz):
            # first layer
            layer = 1 
            if (jj <= horizons[ii,layer]):                    
                velocity_layer[jj,ii] = layers_value[layer-1]                

    # 2nd layer
    for ii in range(0,Nx):
        for jj in range(0,Nz):   
            layer =2
            if (jj > horizons[ii,layer-1] and jj <= horizons[ii,layer]):                    
                velocity_layer[jj,ii] = layers_value[layer-1]

    # 3rd layer
    for ii in range(0,Nx):
        for jj in range(0,Nz):
        # first layer        
            layer =3
            if (jj > horizons[ii,layer-1] and jj <= horizons[ii,layer]):                    
                velocity_layer[jj,ii] = layers_value[layer-1]

    # 4th layer
    for ii in range(0,Nx):
        for jj in range(0,Nz):
        # first layer        
            layer =4
            if (jj > horizons[ii,layer-1] and jj <= horizons[ii,layer]):                    
                velocity_layer[jj,ii] = layers_value[layer-1]


    # 5th layer
    for ii in range(0,Nx):
        for jj in range(0,Nz):
        # first layer        
            layer =5
            if (jj > horizons[ii,layer-1] and jj <= horizons[ii,layer]):                    
                velocity_layer[jj,ii] = layers_value[layer-1]

    # 6th layer
    for ii in range(0,Nx):
        for jj in range(0,Nz):
        # first layer        
            layer =6
            if (jj > horizons[ii,layer-1] and jj <= horizons[ii,layer]):                    
                velocity_layer[jj,ii] = layers_value[layer-1]          
                            

    # last layer
    for ii in range(0,Nx):
        for jj in range(0,Nz):
        # first layer        
            layer =7
            if (jj > horizons[ii,layer-1]):
                velocity_layer[jj,ii] = layers_value[layer-1]

    pl.figure()
    pl.imshow(velocity_layer,cmap='jet')
    

#     # Create velocity model using layers
#     elif (N_horizon >= 2):
#         for ii in range(0,Nx):
#             for jj in range(0,Nz):
#                 # first layer
#                 layer = 0 
#                 if (jj <= horizon[ii,1,0]):
#                     velocity_layer[jj,ii] = layers_min[layer]
#                 # Last layer
#                 layer = N_layers - 1 
#                 if (jj > horizon[ii,1,N_horizon-1]):
#                     velocity_layer[jj,ii] = layers_min[layer]
#                 # second up to N-1 layer    
#                 else:
#                     for layer in range(1,N_layers-1):                    
#                         if ((jj > horizon[ii,1,layer-1]) and jj <= horizon[ii,1,layer]):                        
#                             velocity_layer[jj,ii] = layers_min[layer]


    # # check layered velocity model
    # qc.plotmatrix(velocity_layer,'jet')    
# for layer in range(0,N_layers-1):
#     pl.plot(horizon[:,0,layer],horizon[:,1,layer])

    pl.show()  
    
#     outfile = "xline_buzios_1001z_1202x_layered2.bin"
#     qc.savebinaryfile(Nz,Nx,velocity_layer,outfile)

    

