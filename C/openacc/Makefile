# Makefile to compile and run
# Some tips
# ‘%’ pega o stem (tronco) do nome
# $@ pega o nome do target e 
# $< pega o nome do primeiro pré-requisito
# $^ para listar todos os pré-requisitos do targe
#
# source:https://www.embarcados.com.br/introducao-ao-makefile/

# Compilation parameters
include ../../Makefile.inc
COMPILER=$(CC_GPU)
FLAGS=$(C_GPU_FLAGS)

# Compilation and linking
all: $(C_EXE)

%.exe: %.o
	$(COMPILER) -o $@ $< $(FLAGS)
	
%.o: %.c
	$(COMPILER) -o $@ $< -c $(FLAGS) 

clean:
	@echo "Removing temporaly files ..."
	rm -rf *.o *.mod *~ *.bin *.exe

run:
	@echo "Running program ... "
	./SeismicModeling2D.exe

help:
	@echo $(C_SOURCE) $(OBJ) $(EXE)
	@echo "Makefile"