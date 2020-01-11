# Makefile to compile and run
# Some tips
# ‘%’ pega o stem (tronco) do nome
# $@ pega o nome do target e 
# $< pega o nome do primeiro pré-requisito
# $^ para listar todos os pré-requisitos do targe
#
# source:https://www.embarcados.com.br/introducao-ao-makefile/

# Compilation parameters
include ../Makefile.inc
COMPILER=$(FC)
FLAGS=$(F_FLAGS)

# Compilation and linking
all: $(F_EXE)

%.exe: %.o
	$(COMPILER) -o $@ $< $(FLAGS)
	
%.o: %.f90
	$(COMPILER) -o $@ $< -c $(FLAGS)

clean:
	@echo "Removing temporaly files ..."
	rm -rf *.o *.mod *~ *.bin *.exe

run:
	@echo "Running program ... "
	./SeismicModeling2DVector.exe

help:	
	@echo "Makefile do arquivo modeling_basic.f90"
