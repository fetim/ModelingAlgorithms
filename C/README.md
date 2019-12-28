# Some instructions to show the results

## Using seismic unix package

plot a image:

```ximage n1=numberofsamples < seismogram.bin perc=99 &```

show a movie:

```xmovie n1=dim1 n2=dim2 loop=1 sleep=1 < snapshots.bin &```
