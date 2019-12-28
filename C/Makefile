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
CC=gcc

#Flags for compiler
CC_FLAGS=-c        #\ # flag to create object files
#        -otherflag		  #\ # use backslash to break lines

# ‘%’ pega o stem (tronco) do nome
# $@ pega o nome do target e 
# $< pega o nome do primeiro pré-requisito
# $^ para listar todos os pré-requisitos do targe
#
# Compilation and linking
#
all: $(EXE)
%.exe: %.o
	$(CC) -o $@ $< -lm
	
%.o: %.c
	$(CC) -o $@ $< $(CC_FLAGS) -lm

clean:
	@echo "Removing temporaly files ..."
	rm -rf *.o *.mod *~ *.bin *.exe

run:
	@echo "Running program ... "
	./SeismicModeling2D.exe

help:
	@echo $(C_SOURCE) $(OBJ) $(EXE)
	@echo "Makefile"
