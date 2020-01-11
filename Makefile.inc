# Options for Makefile defined here
# source:https://www.embarcados.com.br/introducao-ao-makefile/

# .C files
C_SRC=$(wildcard *.c)

# .mod files
MOD_SRC=$(wildcard *.mod)

# Object files
OBJ=$(C_SRC:.c=.o)

# Executable files
C_EXE=$(C_SRC:.c=.exe)

# Compiler
CC=gcc
CC_GPU=pgcc
FC=gfortran
FC_GPU=pgfortran

export OMP_NUM_THREADS=11

#Flags for compiler
C_FLAGS=-lm \
		 -fopenmp \
		 -std=c99
C_GPU_FLAGS=-lm \
			-acc -Minfo
F_FLAGS=-O3 \
		 -fopenmp
		 		 
F_GPU_FLAGS=-acc -Minfo

# .f90 files
F_SRC=$(wildcard *.f90)

# Executable files
F_EXE=$(F_SRC:.f90=.exe)


