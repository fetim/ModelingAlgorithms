#!/opt/anaconda3/bin/python
import matplotlib.pyplot as pl
import numpy as np

def readbinaryfile(dim1,dim2,filename):
      """
      readbinaryfile - Functions that read a binary file.
      Usage
      Input:
      dim1     = Number of sample of 1st Dimension
      dim2     = Number of sample of 2nd Dimension
      filename = path of binary file     
      """      
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
      outdata = data.astype('float32')
      outdata.T.reshape(dim1*dim2).tofile(filename)

def plotmatrix(matrix,colormap):
      """
      Plot a 2D matrix read from binary file.
      """      
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

def smooth1D(vector, box_pts):
      box = np.ones(box_pts)/box_pts
      aux = np.ones(box_pts)*vector[0]
      aux = np.append(aux,vector)
      aux = np.append(aux,np.ones(box_pts)*vector[-1])
      aux = np.convolve(aux, box, mode='same')
      vector_smooth = aux[box_pts:-box_pts]

      return vector_smooth

if __name__ =="__main__":
      '''
      Ambiente para teste das funções criadas      
      '''


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