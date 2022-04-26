#------------------------------------------------------------------------------------------------#
#                             MAKEFILE FOR RUNNING COMBUSTION CASE                               #
#------------------------------------------------------------------------------------------------#

# Codes Names
MainCode = Main_Combustion
ExecutableName = Run_Combustion.sh

Class_1 = Memory
Class_2 = ReadData
Class_3 = Parallel
Class_4 = Mesher
Class_5 = PostProcessing
Class_6 = Solver

Directorio_Clases = CppCodes/

default: $(MainCode)

Compilador = mpicxx
Optimisation = -O3

# Libreria PETSc
INCLUDE=-I${PETSC_DIR}/${PETSC_ARCH}/include  -I${PETSC_DIR}/include 
LIB=-Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc

$(MainCode): $(MainCode).cpp
	$(Compilador) -std=c++11 $(INCLUDE) $(MainCode).cpp -o $(ExecutableName)  $(Directorio_Clases)$(Class_1).cpp $(Directorio_Clases)$(Class_2).cpp $(Directorio_Clases)$(Class_3).cpp $(Directorio_Clases)$(Class_4).cpp $(LIB) $(Optimisation) 

clean:
	- $(RM) $(CLEAN) $(ExecutableName)

run:	
	mpirun -np 4 ./$(ExecutableName)
