//------------------------------------------------------------------------------------------------//
//                             MAIN FILE FOR COMBUSTION SIMULATION                                //
//------------------------------------------------------------------------------------------------//

static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <chrono>

// PETSc Library
#include <petsc.h>
#include <petscis.h>
#include <petsctime.h>
#include "petscsys.h"
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "HeaderCodes/Memory.h"
#include "HeaderCodes/ReadData.h"
#include "HeaderCodes/Parallel.h"
#include "HeaderCodes/Mesher.h"
#include "HeaderCodes/PostProcessing.h"
#include "HeaderCodes/Solver.h"

using namespace std;

int main(int argc,char **argv){

PetscInitialize(&argc,&argv,(char*)0, help);

cout << endl;

Memory M1;

ReadData R1(M1);
R1.ReadInputs();

Parallel P1(R1);
P1.RunParallel(M1);

Mesher MESH(M1, R1, P1);
MESH.ExecuteMesher(M1);
P1.Get_MesherInformation(M1, MESH.Ix, MESH.Fx, MESH.NX, MESH.NY, MESH.NZ, MESH.HP, MESH.Halo, MESH.NY_ColumnMP, MESH.NY_ColumnMU, MESH.NY_ColumnMV, MESH.NY_ColumnMW);

MPI_Barrier(MPI_COMM_WORLD);

PostProcessing POST1(M1, R1, MESH, P1);

Solver S1(M1, R1, P1, MESH);

S1.RunSolver(M1, R1, P1, MESH, POST1);

PetscFinalize();


}
