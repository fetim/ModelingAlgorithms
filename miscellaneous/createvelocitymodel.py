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
    
    # load horizons
    fundomar = np.loadtxt("2D_horizons/Fundo_do_Mar_XL3825_2D_davidff.txt")
    pos_sal1 = np.loadtxt("2D_horizons/Pos-Sal_1_XL3825_2D_davidff.txt")
    pos_sal2 = np.loadtxt("2D_horizons/Pos-Sal_2_XL3825_2D_davidff.txt")
    pos_sal3 = np.loadtxt("2D_horizons/Pos-Sal_3_XL3825_2D_davidff.txt")
    toposal  = np.loadtxt("2D_horizons/Topo_Sal_XL3825_2D_davidff.txt")
    basesal  = np.loadtxt("2D_horizons/Base_Sal_XL3825_2D_davidff.txt")

    N_horizon = 6
    N_layers = N_horizon + 1
    horizons = np.zeros([Nx,N_layers])

    # Create horizon table
    horizons[:,0]  = np.arange(0,Nx) # trace index
    horizons[:,1]  = np.rint(fundomar[:-2,4]/dz -1)
    horizons[:,2]  = np.rint(pos_sal1[:-1,4]/dz -1)
    horizons[:,3]  = np.rint(pos_sal2[:,4]/dz -1)
    horizons[:,4]  = np.rint(pos_sal3[:-1,4]/dz -1)
    horizons[:,5]  = np.rint(toposal[:,4]/dz -1)
    horizons[:,6]  = np.rint(basesal[:,4]/dz -1)

    # plot in same window
    fig, ax = pl.subplots()
    ax.imshow(velocity_layered,'jet')        
    for layer in range(1,N_layers):
        # pl.plot(horizon[:,0,layer],horizon[:,1,layer])
        ax.plot(horizons[:,0],horizons[:,layer])     

    pl.colorbar
    pl.show(block=False)  

    # set velocity for each layer        
    horizon_vel = np.zeros([N_layers,1])
    horizon_vel[0] = 1500
    horizon_vel[1] = 1800
    horizon_vel[2] = 2300
    horizon_vel[3] = 2800
    horizon_vel[4] = 3700
    horizon_vel[5] = 4500
    horizon_vel[6] = 5500

    ##### define velocity P
    velocity_horizon = np.zeros([Nz,Nx])
    # Edit first layer of the model
    print("filling layer", layer)
    for ii in range(0,Nx):
        for jj in range(0,Nz):
            # first layer
            layer = 1 
            if (jj <= horizons[ii,layer]):                    
                velocity_horizon[jj,ii] = horizon_vel[layer-1]                

    # 2nd layer - (N-1)th layer
    for layer in range(2,N_layers):
        print("filling layer", layer)
        for ii in range(0,Nx):
            for jj in range(0,Nz):   
                #layer =2
                if (jj > horizons[ii,layer-1] and jj <= horizons[ii,layer]):                    
                    velocity_horizon[jj,ii] = horizon_vel[layer-1]  

    # last layer
    layer = N_layers
    print("filling layer", layer)
    for ii in range(0,Nx):
        for jj in range(0,Nz):
            if (jj > horizons[ii,layer-1]):
                velocity_horizon[jj,ii] = horizon_vel[layer-1]

    # plot in same window
    fig, ax = pl.subplots()
    im = ax.imshow(velocity_horizon,cmap='jet')

    for layer in range(1,N_layers):
            # pl.plot(horizon[:,0,layer],horizon[:,1,layer])
            ax.plot(horizons[:,0],horizons[:,layer])     

    fig.colorbar(im,ax=ax)
    pl.show(block=False) 

    ##### define Vs
    velocity_horizon_vs = velocity_horizon*0.5

    # Edit water layer for vs
    print("filling layer", layer)
    for ii in range(0,Nx):
        for jj in range(0,Nz):
            # first layer
            layer = 1 
            if (jj <= horizons[ii,layer]):                    
                velocity_horizon_vs[jj,ii] = 0

    # plot in same window
    fig, ax = pl.subplots()
    im = ax.imshow(velocity_horizon_vs,cmap='jet')

    for layer in range(1,N_layers):                
            ax.plot(horizons[:,0],horizons[:,layer])     

    fig.colorbar(im,ax=ax)
    pl.show(block=False) 

    # define rho (Gardner)                
    velocity_horizon_rho = 310*np.power(velocity_horizon,0.25)

    # Edit water layer for rho
    print("filling layer", layer)
    for ii in range(0,Nx):
        for jj in range(0,Nz):
            # first layer
            layer = 1 
            if (jj <= horizons[ii,layer]):                    
                velocity_horizon_rho[jj,ii] = 1000

    fig, ax = pl.subplots()
    im = ax.imshow(velocity_horizon_rho,cmap='jet')

    for layer in range(1,N_layers):            
            ax.plot(horizons[:,0],horizons[:,layer])     

    fig.colorbar(im,ax=ax)
    pl.show(block=False) 


    # write velocity P
    outfile = name + '_layer_cake_vp.bin'
    qc.savebinaryfile(Nz,Nx,velocity_horizon,outfile)

    # write velocity S
    outfile = name + '_layer_cake_vs.bin'
    qc.savebinaryfile(Nz,Nx,velocity_horizon_vs,outfile)

    # write density
    outfile = name + '_layer_cake_rho.bin'
    qc.savebinaryfile(Nz,Nx,velocity_horizon_rho,outfile)

    

