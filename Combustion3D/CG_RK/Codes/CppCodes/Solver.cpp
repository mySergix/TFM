//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR CFD SOLVER CLASS                                      //
//------------------------------------------------------------------------------------------------//

#include <chrono>

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"
#include "../HeaderCodes/PostProcessing.h"
#include "../HeaderCodes/Solver.h"

#define N_Species 5
#define N_Reactions 1

Solver::Solver(Memory M1, ReadData R1, Parallel P1, Mesher MESH){

    //Datos Numéricos del problema
	Problema = MESH.Problema;

	NX = MESH.NX;
	NY = MESH.NY;
	NZ = MESH.NZ;

	NX_1 = MESH.NX_1;
	NX_2 = MESH.NX_2;
	NX_3 = MESH.NX_3;

	NY_1 = MESH.NY_1;
	NY_2 = MESH.NY_2;
	NY_3 = MESH.NY_3;
	NY_4 = MESH.NY_4;
	
	Halo = 2;
	HP = 2;

	ProcessNodesP = MESH.ProcessNodesP;
	TotalNodesP = MESH.TotalNodesP;

	// Geometry Data
    Width_Inlet = MESH.Width_Inlet; 
	Width_Slit = MESH.Width_Slit;
	Burner_Wall = MESH.Burner_Wall;
	Symmetry_Burner = MESH.Symmetry_Burner;
	Height_Slit = MESH.Height_Slit;
	Height_Burner =  MESH.Height_Burner;
	Width_Burner = MESH.Width_Burner;

	FlameFront = MESH.FlameFront;

	X_1 = Width_Inlet;
	X_2 = Width_Slit;
	X_3 = Burner_Wall;

	Y_1 = Symmetry_Burner - Height_Burner;
	Y_2 = Height_Slit;
	Y_3 = FlameFront;
	Y_4 = Height_Burner - Height_Slit - FlameFront;

	Xdomain = X_1 + X_2 + X_3;
	if (Problema == "Premixed"){ Ydomain = Y_1 + Y_2 + Y_3 + Y_4; }
	else if (Problema == "NonPremixed"){ Ydomain = Y_2 + Y_3 + Y_4; }
	Zdomain = 0.002;

	//Datos necesarios para computación paralela
	Rango = P1.Rango;
	Procesos = P1.Procesos;
	Ix = M1.AllocateInt(Procesos, 1, 1, 1);
    Fx = M1.AllocateInt(Procesos, 1, 1, 1);

    for (int i = 0; i < Procesos; i++){
        Ix[i] = P1.Ix[i];
        Fx[i] = P1.Fx[i];
    }

	Int_Left = false;
    Int_Right = false;

	CourantFactor = 0.80;
	// Courant Factor -> Driven ----- RK3 ->(0.80 max) ----- RK4 ->(1.10 max) 
	//				  -> Differentially ----- RK3 ->( max) ----- RK4 ->(0.70 max) 

    // Datos Físicos del Problema
	Rho = 1.0;
	Uref = 1.0;
	Reynolds = 10;
	
	U.Gravity = 0.0;
	V.Gravity = 0.0;
	W.Gravity = 0.0;

	Twalls_IntLeft = 298.0; // K
    Twalls_IntRight = 298.0; // K
    Twalls_Slit = 298.0; // K
    Twalls_Burner = 298.0; // K

    T_FlowInlet = 298.0; // K

    ConvergenciaGS = R1.ProblemData[3];
	ConvergenciaGlobal = R1.ProblemData[4];

	R_ideal = 8.31447; // J/mol K

	printf("Process number %d of %d with Ix = %d and Fx = %d, LocalNodes: %d and GlobalNodes: %d \n", Rango, Procesos, Ix[Rango], Fx[Rango], ProcessNodesP, TotalNodesP);
	
}

// Files of the class (CFD)
#include "Matrix_Index.cpp"
#include "Solver_CFD_Memory.cpp"
#include "Solver_CFD_Utilities.cpp"
#include "Solver_CFD_FlowProperties.cpp"
#include "Solver_CFD_BoundaryConditions_Premixed.cpp"
#include "Solver_CFD_BoundaryConditions_Premixed_T.cpp"
//#include "Solver_CFD_BoundaryConditions_NonPremixed.cpp"
//#include "Solver_CFD_BoundaryConditions_NonPremixed_T.cpp"
#include "Solver_CFD_PoissonCoeffs.cpp"
#include "Solver_CFD_MomentumDiffusion.cpp"
#include "Solver_CFD_MomentumConvection.cpp"
#include "Solver_CFD_Energy.cpp"
//#include "Solver_CFD_RK3_Integration.cpp"
#include "Solver_CFD_RK4_Integration.cpp"

#include "Solver_Species_Utilities.cpp"
#include "Solver_Species_Diffusion.cpp"
#include "Solver_Species_Convection.cpp"
#include "Solver_Species_Memory.cpp"
#include "Solver_Species_Read.cpp"
#include "Solver_Species_JANAF.cpp"
#include "Solver_Species_BoundaryConditions_Premixed.cpp"
//#include "Solver_Species_BoundaryConditions_NonPremixed.cpp"
#include "Solver_RunSolver.cpp"

