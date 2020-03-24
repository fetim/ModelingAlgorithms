# ModelingAlgorithms
Repository with acoustic wave equation in some programming languages.

The time of processing using openacc was reduced by 98% using the C algorithm to solve the wave equation. The comparison was made with the serial algorithm.

Instructions of use:

execute to compile:
```
make  
```

in the parameters folder set the file ```2D_acoustic_modeling.dat``` :

```
300    # Nx - Number of horizontal grid points
200    # Nz - Number of vertical grid points
10     # dx - horizontal sample rate
10     # dz - vertical sample rate         
2001   # Nt - Number of time samples
1.0e-3 # dt - time sample rate      
3      # Nshot - Number of shot   
30     #     
```

Choose the programing language and enter in the folder. For example:

``` cd C/serial ```

Execute the modeling algorithm with:

``` make run ```

Check the results with ```ximage``` and ```xmovie```

Seismogram:
``` ximage n1=Nt < seismogram.bin perc=99 &```

Snapshots:
``` xmovie n1=Nz n2=Nx loop=1 sleep=1 < snapshots.bin & ```

