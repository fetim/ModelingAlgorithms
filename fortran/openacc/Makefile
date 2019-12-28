# source:https://www.embarcados.com.br/introducao-ao-makefile/
# author: @fetim
# Name of project
PROJ_NAME=Basicmodeling

# .f90 files
F90_SOURCE=$(wildcard *.f90)

# .mod files
MOD_SOURCE=$(wildcard *.mod)

# Object files
OBJ=$(F90_SOURCE:.f90=.o)

# Executable files
EXE=$(F90_SOURCE:.f90=.exe)

# Compiler
FC=pgfortran

#Flags for compiler
FC_FLAGS=-c        #\ # flag to create object files
#        -otherflag		  #\ # use backslash to break lines

PARALLEL_FLAG=-acc -Minfo
# ‘%’ pega o stem (tronco) do nome
# $@ pega o nome do target e 
# $< pega o nome do primeiro pré-requisito
# $^ para listar todos os pré-requisitos do targe
#
# Compilation and linking
#
all: $(EXE)
%.exe: %.o
	$(FC) -o $@ $< $(PARALLEL_FLAG)
	
%.o: %.f90
	$(FC) -o $@ $< $(FC_FLAGS) $(PARALLEL_FLAG)

clean:
	@echo "Removing temporaly files ..."
	rm -rf *.o *.mod *~ *.bin *.exe

run:
	@echo "Running program ... "
	./SeismicModeling2DVector.exe

help:
	@echo $(F90_SOURCE) $(OBJ) $(EXE)
	@echo "Makefile do arquivo modeling_basic.f90"
