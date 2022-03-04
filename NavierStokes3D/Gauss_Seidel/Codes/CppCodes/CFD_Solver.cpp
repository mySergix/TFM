//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR CFD SOLVER CLASS                                      //
//------------------------------------------------------------------------------------------------//

#include <chrono>

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"
#include "../HeaderCodes/PostProcessing.h"
#include "../HeaderCodes/CFD_Solver.h"

CFD_Solver::CFD_Solver(Memory M1, ReadData R1, Parallel P1){

    //Datos Numéricos del problema
	Problema = R1.ProblemNumericalData[0];

	NX = R1.ProblemNumericalData[2];
	NY = R1.ProblemNumericalData[3];
	NZ = R1.ProblemNumericalData[4];
	Halo = 2;
	HP = 2;

    //Datos Geométricos del problma
	Xdominio = R1.GeometryData[0];
	Ydominio = R1.GeometryData[1];
	Zdominio = R1.GeometryData[2];

	//Datos necesarios para computación paralela
	Rango = P1.Rango;
	Procesos = P1.Procesos;
	Ix = M1.AllocateInt(Procesos, 1, 1, 1);
    Fx = M1.AllocateInt(Procesos, 1, 1, 1);

    for (int i = 0; i < Procesos; i++){
        Ix[i] = P1.Ix[i];
        Fx[i] = P1.Fx[i];
    }

    //Datos Físicos del Problema
	Rho = R1.ProblemPhysicalData[0];
	Uref = R1.ProblemPhysicalData[1];
	Reynolds = R1.ProblemPhysicalData[2];
	
	Rayleigh = R1.ProblemPhysicalData[3];
	Cp = R1.ProblemPhysicalData[4];
	Prandtl = R1.ProblemPhysicalData[5];

	U.Gravity = R1.ProblemPhysicalData[6];
	V.Gravity = R1.ProblemPhysicalData[7];
	W.Gravity = R1.ProblemPhysicalData[8];

	Tleft = R1.ProblemPhysicalData[9];
	Tright = R1.ProblemPhysicalData[10];

	if (Problema == 1){
		mu = (Uref * Xdominio) / Reynolds;
	}
	else if (Problema == 2){
		To = abs(Tleft + Tright) / 2.0; 
		Beta = 1.0 / To;	
		mu = sqrt( (pow(Rho,2.0) * abs(V.Gravity) * pow(Xdominio,3) * Beta * abs(Tleft - Tright) * Prandtl) / Rayleigh);
		K = (Cp * mu) / Prandtl;
	}
	
    ConvergenciaGS = R1.ProblemData[3];
	ConvergenciaGlobal = R1.ProblemData[4];

}

// Files of the class (CFD)
#include "Matrix_Index.cpp"
#include "CFD_Solver_Memory.cpp"
#include "CFD_Solver_Utilities.cpp"
#include "CFD_Solver_BoundaryConditions_Driven.cpp"
//#include "CFD_Solver_BoundaryConditions_Differentially.cpp"
#include "CFD_Solver_PoissonCoeffs.cpp"
#include "CFD_Solver_MomentumDiffusion.cpp"
#include "CFD_Solver_MomentumConvection.cpp"
#include "CFD_Solver_Energy.cpp"
#include "CFD_Solver_RunSolver.cpp"