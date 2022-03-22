#------------------------------------------------------------------------------------------------#
#                             MAKEFILE FOR RUNNING MASS TRANSPORT CASE                           #
#------------------------------------------------------------------------------------------------#

# Codes Names
MainCode = Main_MassTransport
ExecutableName = Run_MassTransport

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
ChronoLib_Flag = -std=c++11

# Libreria PETSc
INCLUDE=-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include 
LIB=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc

$(MainCode): $(MainCode).cpp
	$(Compilador) -g -std=c++11 -o $(ExecutableName) $(MainCode).cpp $(Directorio_Clases)$(Class_1).cpp $(Directorio_Clases)$(Class_2).cpp $(Directorio_Clases)$(Class_3).cpp $(Directorio_Clases)$(Class_4).cpp $(Directorio_Clases)$(Class_5).cpp $(Directorio_Clases)$(Class_6).cpp $(INCLUDE) $(LIB) $(Optimisation)

clean:
	- $(RM) $(CLEAN) $(ExecutableName)

run:	
	mpirun -np 4 ./$(ExecutableName)