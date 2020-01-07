# source:https://www.embarcados.com.br/introducao-ao-makefile/
# author: @fetim
# Name of project
PROJ_NAME=Basicmodeling

# .C files
C_SOURCE=$(wildcard *.c)

# .mod files
MOD_SOURCE=$(wildcard *.mod)

# Object files
OBJ=$(C_SOURCE:.c=.o)

# Executable files
EXE=$(C_SOURCE:.c=.exe)

# Compiler
CC=pgcc

#Flags for compiler
CC_FLAGS=-c       #\ # flag to create object files
#        -otherflag		  #\ # use backslash to break lines

# ‘%’ pega o stem (tronco) do nome
# $@ pega o nome do target e 
# $< pega o nome do primeiro pré-requisito
# $^ para listar todos os pré-requisitos do targe
#
# Compilation and linking
#
all: $(EXE)
# $(EXE): $(OBJ)
# 	@echo ""
# 	@echo "main program generation ..."
# 	$(CC) -o $@ $< $(PARALLEL_FLAG)

%.exe: %.o
	$(CC) -o $@ $< -lm -fopenmp
	
%.o: %.c
#	 @echo ""
#	 @echo "object files generation ..."
	$(CC) -o $@ $< $(CC_FLAGS)-lm -fopenmp

clean:
	@echo "Removing temporaly files ..."
	rm -rf *.o *.mod *~ *.bin *.exe

run:
	@echo "Running program ... "
	./SeismicModeling2D.exe

help:
	@echo $(C_SOURCE) $(OBJ) $(EXE)
	@echo "Makefile"
