//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR CFD SOLVER CLASS                                   //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <chrono>
#include <algorithm>

#include "petsc.h"

#define N_Species 5
#define N_Reactions 4

using namespace std;

// Forward declaration of classes

class Solver{
    private:

    public:

        Solver(Memory, ReadData, Parallel, Mesher);

        // Problem type
		string Problema;

		int NX;
		int NY; 
		int NZ;

		int NX_1, NX_2, NX_3;
		int NY_1, NY_2, NY_3, NY_4;

        int ProcessNodesP;
        int TotalNodesP;
        
        // Geometry Data
		double Width_Inlet;
		double Width_Slit;
		double Burner_Wall;
		double Symmetry_Burner;
		double Height_Slit;
		double Height_Burner;
		double Width_Burner;

		double FlameFront;
		
		double X_1, X_2, X_3;
		double Y_1, Y_2, Y_3, Y_4;

		double Xdomain, Ydomain, Zdomain;

        string EsquemaLargo;
		string EsquemaCorto;

        // Parameters for parallel computing
		int Rango;
		int Procesos;
        int Halo;
        int HP;

		int *Ix;
		int *Fx;
		
        bool Int_Left;
        bool Int_Right;
        
        // Physical Data
		double Rho;
		double Uref;
        double Uref_NonPremixed;
		double Reynolds;

        double Twalls_IntLeft;
        double Twalls_IntRight;
        double Twalls_Slit;
        double Twalls_Burner;

        double T_FlowInlet;
        
		double Rayleigh;
		double Cp;
		double Prandtl;
        double K;
        
        double CourantFactor;
        double DiffusiveDeltaT;
		double DeltaT;
        double nu;
        double mu;
    
        double ConvergenciaGS;
		double MaxDiffGS;

		double ConvergenciaGlobal;
		double MaxDiffGlobal;
        
        // PETSc Data
        PetscInt    NNZ;
        PetscInt	*Col_Ind;
	    PetscInt	*Row_Ptr;
	    PetscScalar *Val_Laplacian;

        PetscInt    *RHS_Ind;
        PetscScalar *RHS;

        PetscScalar *X_Sol_Array;
        
        Vec         B_RHS, X_Sol;
        Mat         A_Matrix;

        // Flow Properties
        double *mu_visc;
        double *K_Thermal;
        double *Cp_Heat;

        string *SpeciesString;
        double R_ideal;

        struct Velocity_Struct
        {
            double *Predictor;

            double *Pres;
            double *Fut;

            double *Convective;
            double *Diffusive;

            double *Bottom;
            double *Top;

            double *Here;
            double *There;

            double *Left;
            double *Right;

            double *I_Left; // Internal Left
            double *I_Right; // Internal Right

            double *Boussinesq;

            double *K1;
            double *K2;
            double *K3;
            double *K4;

            double *New_Velocity;

            double Gravity;
        };

        struct Energy_Struct
        {
            double *Pres;
            double *Fut;

            double *Convective;
            double *Diffusive;
            double *Reactive;

            double *ContributionPast;
            double *ContributionPres;

            double *Bottom;
            double *Top;

            double *Here;
            double *There;

            double *Left;
            double *Right;

            double *I_Left; // Internal Left
            double *I_Right; // Internal Right
        }; 

        struct Poisson_Coeffs
        {
            double *aw;
            double *ae;

            double *as;
            double *an;

            double *ah;
            double *at;

            double *ap;
        };

        struct Pressure_Struct
        {
            double *Pres;
            double *Sup;
        };

        struct Global_Struct
        {
            double *P;
            double *U;
            double *V;
            double *W;

            double *T;
        };

        struct Runge_Kutta
        {
            double c1, c2, c3, c4;
            
            double a11, a12, a13, a14;
            double a21, a22, a23, a24;
            double a31, a32, a33, a34;
            double a41, a42, a43, a44;

            double b1, b2, b3, b4;
        };

        struct Species_Struct
        {
            string Name; // Name of the species
            double Wmolar; // Molar weight of the species
            double Epsilon; // Characteristic Lennard-Jones energy
            double Schmidt; // Schmidt number of the specie
            double Lewis; // Lewis number of the specie

            double sigma;

            double InitialY; // Initial Mass Fraction in the domain

            double *Y_Pres;
            double *Y_Fut;

            double *Convective;
            double *Diffusive;

            double *ContributionPast;
            double *ContributionPres;

            double *D_ab;

            double *Bottom;
            double *Top;

            double *Here;
            double *There;

            double *Left;
            double *Right;

            // JANAF Terms
            double *Cp_coeff;
            double *h_coeff;
            double *mu_coeff;
            double *lambda_coeff;

            // Chemical Terms
            double *wk;

            double *Global;
        };

        struct Reactions_Struct
        {
            double *Kf;
            double *Kr;

            double *A;
            double *Beta;
            double *EA;

            double *DeltaS;
            double *DeltaH;

        };

        struct Velocity_Struct U;
        struct Velocity_Struct V;
        struct Velocity_Struct W;

        struct Pressure_Struct P;

        struct Energy_Struct T;

        struct Poisson_Coeffs A;
        
        struct Global_Struct Global;

        struct Runge_Kutta RK;

        struct Species_Struct Species[N_Species];

        struct Reactions_Struct Reactions;

        // Class Functions

            // Memory Allocation
            void Allocate_PoissonCoeffsMemory(Memory);
            void Allocate_VelocitiesPartMemory(Memory, Velocity_Struct&, int, int, int);
            void Allocate_VelocitiesBoundaryConditionsMemory(Memory, Velocity_Struct&, int, int, int);
            void Allocate_PressureMemory(Memory);
            void Allocate_VelocitiesMemory(Memory);
            void Allocate_EnergyMemory(Memory);
            void Allocate_GlobalMemory(Memory);

            // Memory Delete
            void Delete_EnergyMemory(Energy_Struct&);
            void Delete_VelocityMemory(Velocity_Struct&);
            void Delete_PoissonMemory(Poisson_Coeffs&);
            void Delete_SolverMemory();

            // Utilities
            void Get_InitialConditions(Mesher);
            void Get_ProcessesInternalCheck();
            void Get_DiffusiveTimeStep(Mesher);
            void Get_StepTime(Mesher);
            inline double ConvectiveScheme(double, double, double, double, double, double, double, double, double, double, string);  
            void Get_PredictorsDivergence(Mesher);
            void Get_CorrectedVelocities(Mesher, Parallel, double*, double*, double*, double);
            void Get_Velocities(Mesher, Parallel);
            void Get_PressureHalos();
            void Get_Stop(Mesher);
            void Get_Update(Mesher);

            // Flow Properties
            void Get_DynamicViscosity(Mesher);
            void Get_ThermalConductivity(Mesher);
            void Get_CpHeat(Mesher);

            // Runge Kutta Temporal Integration
            void Get_RK_Coefficients(Runge_Kutta&);
            void Get_RK_VelocityContributions(Mesher, double*, double*, double*);
            void Get_RK_Integration(Mesher, Parallel);

            // Boundary Conditions
            void Get_StaticBoundaryConditions_Velocities(Mesher);
            void Get_UpdateBoundaryConditions_Velocities(Mesher, double*, double*, double*);
            void Get_UpdateBoundaryConditions_PredictorVelocities(Mesher);
            void Get_StaticHalos_Velocity(Mesher, double*, double*, double*);
            void Get_UpdateHalos_Velocity(Mesher, double*, double*, double*);
            
            void Get_StaticBoundaryConditions_Temperatures(Mesher);
            void Get_UpdateBoundaryConditions_Temperatures(Mesher, double*);
            void Get_StaticHalos_Temperatures(Mesher, double*);
            void Get_UpdateHalos_Temperatures(Mesher, double*);

            // Poisson Coefficients
            void Get_PoissonCoefficients(Mesher);
            void Get_NonZero_NumberElements();
            int Get_GlobalIndCoefficient(Mesher, int, int, int);
            int Get_LocalIndCoefficient(Mesher, int, int, int);
            void Get_CSR_LaplacianMatrix(Mesher);
            void Get_RHS_VectorIndex(Mesher);
            void Get_LocalSolution(Mesher);

            // Momentum Diffusion
            void Get_DiffusionU(Mesher, double*);
            void Get_DiffusionV(Mesher, double*);
            void Get_DiffusionW(Mesher, double*);

            // Momentum Convection
            void Get_ConvectionU(Mesher, double*, double*, double*);
            void Get_ConvectionV(Mesher, double*, double*, double*);
            void Get_ConvectionW(Mesher, double*, double*, double*);

            // Energy Equation
            void Get_DiffusionEnergy(Mesher);
            void Get_ConvectionEnergy(Mesher);
            void Get_EnergyContributions(Mesher);
            void Get_NewTemperatures(Mesher);

            // Species Transport Equation
            void Allocate_StructSpecies(Memory);
            void Delete_StructSpecies();

            void Get_Species_StaticBoundaryConditions(Mesher);
            void Get_Species_UpdateBoundaryConditions(Mesher);

            void Get_Species_StaticHalos();
            void Get_Species_UpdateHalos(Mesher);

            void Get_Species_Convection(Mesher);

            void Get_Species_DiffusionCoefficients(Mesher);
            void Get_Species_Diffusion(Mesher);

            void Get_Species_InitialConditions(Mesher);
            void Get_SpeciesDiffusion_TimeStep(Mesher);
            void Get_SpeciesContributions(Mesher);
            void Get_Species_MassFraction(Mesher);
            void Get_Species_MassConservation(Mesher);
            void Get_Species_Update(Mesher);

            // Combustion Data Input
            void Read_SpeciesName(string);
            void Read_AllSpeciesData();
            void Read_SpeciesInformation(Species_Struct&, string);

            // Combustion JANAF Calculations
            double JANAF_CpSpecie(double, int);
            double JANAF_LambdaSpecies(double, int);
            double JANAF_CpHeat(double, int, int, int);
            double JANAF_AbsEnthalpy_Specie(int, double);
            double JANAF_AbsEnthalpy_Specie_Mix(double, int, int, int);
            double JANAF_DynViscosity(double, int, int, int);
            double JANAF_ThermalCond(double, int, int, int);

            void Get_DynamicViscosity();
            void Get_ThermalConductivity();
            void Get_CpHeat();

            // Combustion-Related Calculations
            void Allocate_StructReactions(Memory);
            void Get_ReactionsEnergy(Mesher);

            void Get_Reactions_OneStepCH4(Mesher);
            void Get_Reactions_FourStepCH4(Mesher);

            // Run Solver
            void RunSolver(Memory, ReadData, Parallel, Mesher, PostProcessing);
            
};