//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR CFD SOLVER CLASS                                      //
//------------------------------------------------------------------------------------------------//

#include <chrono>

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"
//#include "../HeaderCodes/PostProcessing.h"
#include "../HeaderCodes/Solver.h"

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

	CourantFactor = 0.70;
	// Courant Factor -> Driven ----- RK3 ->(0.80 max) ----- RK4 ->(1.10 max) 
	//				  -> Differentially ----- RK3 ->( max) ----- RK4 ->(0.70 max) 

    // Datos Físicos del Problema
	Rho = ;
	Uref = ;
	Reynolds = ;
	
	Prandtl = R1.ProblemPhysicalData[5];

	U.Gravity = ;
	V.Gravity = ;
	W.Gravity = ;

	mu = (Uref * Xdomain) / Reynolds;
	
    ConvergenciaGS = R1.ProblemData[3];
	ConvergenciaGlobal = R1.ProblemData[4];

}

// Files of the class (CFD)
#include "Matrix_Index.cpp"
#include "Solver_CFD_Memory.cpp"
#include "Solver_CFD_Utilities.cpp"
//#include "Solver_CFD_BoundaryConditions_Driven.cpp"
#include "Solver_CFD_BoundaryConditions_Differentially.cpp"
#include "Solver_CFD_BoundaryConditions_Differentially_T.cpp"
#include "Solver_CFD_PoissonCoeffs.cpp"
#include "Solver_CFD_MomentumDiffusion.cpp"
#include "Solver_CFD_MomentumConvection.cpp"
#include "Solver_CFD_Energy.cpp"
//#include "Solver_CFD_RK3_Integration.cpp"
#include "Solver_CFD_RK4_Integration.cpp"
#include "Solver_RunSolver.cpp"

