//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR CFD SOLVER CLASS                                   //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <chrono>

#include "petsc.h"

#define N_Species 3
#define N_Reactions 4

using namespace std;

// Forward declaration of classes

class Solver{
    private:

    public:

        Solver(Memory, ReadData, Parallel);

        // Problem type
		int Problema;

        // Geometry Data
		double Xdominio;
		double Ydominio;
		double Zdominio;

        // Problem Data
		int NX;
		int NY; 
		int NZ;

        string EsquemaLargo;
		string EsquemaCorto;

        // Parameters for parallel computing
		int Rango;
		int Procesos;
        int Halo;
        int HP;

		int *Ix;
		int *Fx;
		
        // Physical Data
		double Rho;
		double Uref;
		double Reynolds;

		double Rayleigh;
		double Cp;
		double Prandtl;

		double gx;
		double gy;
		double gz;

        double CourantFactor;
        double DiffusiveDeltaT;
		double DeltaT;

        double Tleft;
        double Tright;
    
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

        string *SpeciesString;

        double R_ideal;

        // Channel Mass Transport Study Case Variables
        double Dh;
        double Gamma_Geometry;

        double To;
        double Twater;
		double Beta;
        double Beta_Y;
		double K;
        double mu;

        double Schmidt;
        double MW_Air;
        double MW_Water;

        double w_av;
        double Co;
        double hfg;
        double D_AB;

        double Qs; // Mass Flow

        double *mu_visc;
        double *K_Thermal;
        double *Cp_Heat;

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

            double *Boussinesq;
            double *Boussinesq_Y;

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

            double *Tbottom;
            double *Ybottom;
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

            double *Y_Pres;
            double *Y_Fut;

            double *X; // Molar fraction

            double *Convective;
            double *Diffusive;

            double *ContributionPast;
            double *ContributionPres;

            double *D_AlphaMix;

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

        struct Species_Struct Species[N_Species];

        struct Reactions_Struct Reactions;

        struct Velocity_Struct U;
        struct Velocity_Struct V;
        struct Velocity_Struct W;

        struct Pressure_Struct P;

        struct Energy_Struct T;
        struct Poisson_Coeffs A;
        
        struct Global_Struct Global;

        struct Runge_Kutta RK;


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
            void Get_MassFlowValue(Mesher);
            void Get_DiffusiveTimeStep(Mesher);
            void Get_StepTime(Mesher);
            inline double ConvectiveScheme(double, double, double, double, double, double, double, double, double, double, string);  
            void Get_PredictorsDivergence(Mesher);
            void Get_CorrectedVelocities(Mesher, Parallel, double*, double*, double*);
            void Get_Velocities(Mesher, Parallel);
            void Get_Stop();
            void Get_Update();

            // Runge Kutta Temporal Integration
            void Get_RK_Coefficients(Runge_Kutta&);
            void Get_RK_VelocityContributions(double*, double*, double*);
            void Get_RK_Integration(Mesher, Parallel);

            // Boundary Conditions
            void Get_StaticBoundaryConditions_Velocities(Mesher);
            void Get_UpdateBoundaryConditions_Velocities(double*, double*, double*);
            void Get_UpdateBoundaryConditions_PredictorVelocities();
            void Get_StaticHalos_Velocity(double*, double*, double*);
            void Get_UpdateHalos_Velocity(double*, double*, double*);
            
            void Get_StaticBoundaryConditions_Temperatures(Mesher);
            void Get_UpdateBoundaryConditions_Temperatures(Mesher, double*);
            void Get_StaticHalos_Temperatures(double*);
            void Get_UpdateHalos_Temperatures(double*);

            // Poisson Coefficients
            void Get_PoissonCoefficients(Mesher);
            void Get_NonZero_NumberElements();
            void Get_CSR_LaplacianMatrix();
            void Get_RHS_VectorIndex();
            void Get_LocalSolution();

            // Momentum Diffusion
            void Get_DiffusionU(Mesher, double*);
            void Get_DiffusionV(Mesher, double*);
            void Get_DiffusionW(Mesher, double*);

            // Momentum Convection
            void Get_ConvectionU(Mesher, double*, double*, double*);
            void Get_ConvectionV(Mesher, double*, double*, double*);
            void Get_ConvectionW(Mesher, double*, double*, double*);

            // Energy Equation
            void Get_DiffusionEnergy(Mesher, double*);
            void Get_ConvectionEnergy(Mesher, double*, double*, double*, double*);
            void Get_BoussinesqV(Mesher, double*, double*);
            void Get_EnergyContributions();
            void Get_NewTemperatures(Mesher);

            // Species Transport Equation
            void Allocate_StructSpecies(Memory);
            void Delete_StructSpecies();

            void Get_Species_StaticBoundaryConditions(Mesher);
            void Get_Species_UpdateBoundaryConditions(Mesher);
            void Get_Species_StaticHalos();
            void Get_Species_UpdateHalos();

            void Get_Species_Convection(Mesher);
            void Get_Species_DiffusionCoefficients();
            void Get_Species_Diffusion(Mesher);

            void Get_Species_InitialConditions();
            void Get_Boussinesq_MassFraction(Mesher, double*);
            void Get_SpeciesContributions();
            void Get_Species_MassFraction();
            void Get_Species_MassConservation();
            
            // Combustion Data Input
            void Read_SpeciesName(string);
            void Read_AllSpeciesData();
            void Read_SpeciesInformation(Species_Struct&, string);

            // Combustion JANAF Calculations
            double JANAF_CpSpecie(double, int);
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
            void Get_ReactionsEnergy();

            // Diffusion models
            double Get_Diffusion_Schmidt(int);
            double Get_Diffusion_Lewis(int, double);
            double Get_Diffusion_ChampanEnskog(int, double, double, int, int, int);
            double Get_Diffusion_WilkeLee(int, double, double, int, int, int);

            // Run Solver
            void RunSolver(Memory, ReadData, Parallel, Mesher, PostProcessing);
            
};