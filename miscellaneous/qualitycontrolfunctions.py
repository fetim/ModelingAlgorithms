#!/usr/bin/python
#  Plot 2D binary file.          
#  This script import a 2D binary file
#  and plot as image
#
#  INPUT:  
#  filename      = path with the 2D binary file.
#  dim1          = Number of elements in first dimension
#  dim2          = Number of elements in second dimension
# 
#  OUTPUT: 
#  file          = None;
#  
#  Code Written by Felipe Timoteo and Cintia Queiroz
#                  Last update: 10 May, 2019
# 
#  Copyright (C) 2019 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
#                     Departamento de Geologia e Geofísica
#                     Universidade Federal Fluminense
###############################################################################

def readbinaryfile(dim1,dim2,filename):
      """
      readbinaryfile - Functions that read a binary file.
      Usage
      Input:
      dim1     = Number of sample of 1st Dimension
      dim2     = Number of sample of 2nd Dimension
      filename = path of binary file     
      """      
      import numpy as np
      with open(filename, 'rb') as f:    
            data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
            matrix = np.reshape(data, [dim1,dim2], order='F')
      return matrix

def savebinaryfile(dim1,dim2,data,filename):
      """
      savebinaryfile - Functions that read a binary file.
      Usage
      Input:
      dim1     = Number of sample of 1st Dimension
      dim2     = Number of sample of 2nd Dimension
      data     = 2D array
      filename = path of binary file     
      """      
      import numpy as np

      outdata = data.astype('float32')
      outdata.T.reshape(dim1*dim2).tofile(filename)

def plotmatrix(matrix,colormap):
      """
      Plot a 2D matrix read from binary file.
      """      
      import matplotlib.pyplot as pl
      from mpl_toolkits.axes_grid1 import make_axes_locatable

      pl.figure()
      ax = pl.gca()
      im = pl.imshow(matrix,cmap= colormap,aspect='auto')

      # create an axes on the right side of ax. The width of cax will be 5%
      # # of ax and the padding between cax and ax will be fixed at 0.05 inch.
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)

      pl.colorbar(im, cax=cax)
      # pl.draw(block=False) # drawing figure to be plotted later 
      pl.show(block=False) 

if __name__ =="__main__":
      '''
      Ambiente para teste das funções criadas      
      '''
      import matplotlib.pyplot as pl
      import numpy as np

      import numpy as np

      f = open('../parameters/2D_acoustic_modeling.dat', 'r') # 'r' = read
      parameters = np.genfromtxt(f,delimiter='')

      Nx     = int(parameters[0])
      Nz     = int(parameters[1])
      dx     = parameters[2]
      dz     = parameters[3]
      Nt     = int(parameters[4])
      dt     = parameters[5]
      Nshot  = int(parameters[6])
      fcut   = parameters[7]

      print("Parameters: ")
      print("Nx    = ",Nx  )
      print("Nz    = ",Nz  )
      print("dx    = ",dx  )
      print("dz    = ",dz  )
      print("Nt    = ",Nt  )
      print("dt    = ",dt  )
      print("Nshot = ",Nshot)
      print("fcut  = ",fcut)

      Nsnapshot = 10
      snap      = 7

      # inputname  = '../C/snapshots.bin'
      # snapshot_A = readbinaryfile(Nz,Nx*Nsnapshot,inputname)
      # plotmatrix(snapshot_A[:,snap*Nx:(snap+1)*Nx],'gray')
      
      # inputname  = '../fortran/snapshots.bin'
      # snapshot_B = readbinaryfile(Nz,Nx*Nsnapshot,inputname)
      # plotmatrix(snapshot_B[:,snap*Nx:(snap+1)*Nx],'gray')

      inputname  = '../C/openmp/seismogram.bin'
      seismogram_A = readbinaryfile(Nt,Nx*Nshot,inputname)
      plotmatrix(seismogram_A[:,:],'gray')

      inputname  = '../fortran/openmp/seismogram.bin'
      seismogram_B = readbinaryfile(Nt,Nx*Nshot,inputname)
      plotmatrix(seismogram_B[:,:],'gray')

      pl.show()