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

#define N_Species 3

Solver::Solver(Memory M1, ReadData R1, Parallel P1){

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

	CourantFactor = 0.70;
	// Courant Factor -> Driven ----- RK3 ->(0.80 max) ----- RK4 ->(1.10 max) 
	//				  -> Differentially ----- RK3 ->( max) ----- RK4 ->(0.70 max) 

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

    ConvergenciaGS = R1.ProblemData[3];
	ConvergenciaGlobal = R1.ProblemData[4];

	R_ideal = 8.314;
	
	Dh = 4.0 * (Zdominio * Ydominio) / (2.0 * (Zdominio + Ydominio)); // Hidraulic Diameter
    Gamma_Geometry = Ydominio / Zdominio;

    To = 22.40; // Celsius
    Twater = 16.0; // Celsius

	Beta = 0.00343;
    Beta_Y = 0.513;

    MW_Air = 28.97; // g/mol
    MW_Water = 18.0; // g/mol

	Schmidt = 0.60;
	mu = sqrt((abs(V.Gravity) * Beta * abs(Twater - To) * pow(Rho, 2.0) * pow(Dh, 3.0) * Prandtl) / Rayleigh);
	K = (mu * Cp) / Prandtl;
	D_AB = mu / (Rho * Schmidt);

    w_av = (Reynolds * mu) / (Rho * Dh);
    Co = 0.531 * MW_Water / 1000; //(100 * pow(10, 0.66077 + (7.5 * To) / (237.3 + To)) * (MW_Water / 1000)) / (8.314 * (273.15 + To));
    hfg = (40.7e3) / (MW_Water / 1000.0);

	Qs = 0.0;

}

// Files of the class (CFD)
#include "Matrix_Index.cpp"
#include "Solver_CFD_Memory.cpp"
#include "Solver_CFD_Utilities.cpp"
//#include "Solver_CFD_BoundaryConditions_Driven.cpp"
//#include "Solver_CFD_BoundaryConditions_Differentially.cpp"
//#include "Solver_CFD_BoundaryConditions_Differentially_T.cpp"
#include "Solver_CFD_BoundaryConditions_Channel.cpp"
#include "Solver_CFD_BoundaryConditions_Channel_T.cpp"
#include "Solver_CFD_PoissonCoeffs.cpp"
#include "Solver_CFD_MomentumDiffusion.cpp"
#include "Solver_CFD_MomentumConvection.cpp"
#include "Solver_CFD_Energy.cpp"
//#include "Solver_CFD_RK3_Integration.cpp"
#include "Solver_CFD_RK4_Integration.cpp"
#include "Solver_RunSolver.cpp"

// Files of the class (Species)
#include "Solver_Species_BoundaryConditions.cpp"
#include "Solver_Species_Convection.cpp"
#include "Solver_Species_Diffusion.cpp"
#include "Solver_Species_Memory.cpp"
#include "Solver_Species_Utilities.cpp"

// Files of the class (Combustion)
#include "Solver_Species_Read.cpp"
#include "Solver_Species_JANAF.cpp"