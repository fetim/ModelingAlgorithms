Instructions of use:

execute to compile:
```
make  
```

if doesn't work try (needed gcc compiler or similar):

```gcc -o SeismicModeling2D.exe SeismicModeling2D.c -lm```

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

or

```./Seismogram.exe ```

Check the results with ```ximage``` and ```xmovie``` (seismic unix package)

Seismogram:
``` ximage n1=Nt < seismogram.bin perc=99 &```

Snapshots:
``` xmovie n1=Nz n2=Nx loop=1 sleep=1 < snapshots.bin & ```

plot a image:

```ximage n1=numberofsamples < seismogram.bin perc=99 &```

show a movie:

```xmovie n1=dim1 n2=dim2 loop=1 sleep=1 < snapshots.bin &```
