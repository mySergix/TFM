//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR CFD SOLVER CLASS                                   //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#define N_Species 2

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

        double Dh;
        double gamma_geometry;

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

		double DeltaT;
        double nu;
        double mu;

        double Tleft;
        double Tright;
        
        double To;
		double Beta;
		double K;
    
        double ConvergenciaGS;
		double MaxDiffGS;

		double ConvergenciaGlobal;
		double MaxDiffGlobal;
        
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

            double *K1;
            double *K2;
            double *K3;

            double *New_Velocity;

            double Gravity;
        };

        struct Energy_Struct
        {
            double *Pres;
            double *Fut;

            double *Convective;
            double *Diffusive;

            double *K1;
            double *K2;
            double *K3;

            double *New_Temperature;

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
            double *bp;
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
            double c1, c2, c3;
            
            double a11, a12, a13;
            double a21, a22, a23;
            double a31, a32, a33;

            double b1, b2, b3;
        };

        struct Velocity_Struct U;
        struct Velocity_Struct V;
        struct Velocity_Struct W;

        struct Pressure_Struct P;

        struct Energy_Struct T;
        struct Poisson_Coeffs A;
        
        struct Global_Struct Global;

        struct Runge_Kutta RK3;

        struct Species_Struct
        {
            double *C_Pres;
            double *C_Fut;

            double *ContributionPast;
            double *ContributionPres;

            double *Convective;
            double *Diffusive;

            double *D_ab;

            double *Bottom;
            double *Top;

            double *Here;
            double *There;

            double *Left;
            double *Right;

            double *Global;
        };

        Species_Struct Species[N_Species];

        // Class functions

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
            void Get_StepTime(Mesher, Parallel);
            inline double ConvectiveScheme(double, double, double, double, double, double, double, double, double, double, string);
            void Get_ContributionsPredictors();
            void Get_PredictorsDivergence(Mesher);
            void Get_Velocities(Mesher, Parallel);
            void Get_Stop();
            void Get_Update();
            void Get_ConcentrationBuoyancy(Mesher, int);

            // Boundary Conditions
            void Get_StaticBoundaryConditions_Velocities(Mesher);
            void Get_UpdateBoundaryConditions_Velocities(double*, double*, double*);
            void Get_UpdateBoundaryConditions_PredictorVelocities();
            void Get_StaticHalos_Velocity();
            void Get_UpdateHalos_Velocity(double*, double*, double*);
            
            void Get_StaticBoundaryConditions_Temperatures(Mesher);
            void Get_UpdateBoundaryConditions_Temperatures(double*);
            void Get_StaticHalos_Temperatures();
            void Get_UpdateHalos_Temperatures(double*);

            // Poisson Coefficients
            void Get_PoissonCoefficients(Mesher);
            void Get_GaussSeidel(Parallel);
            void PrintTxt();

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

            // Runge Kutta Temporal Integration
            void Get_RK_Coefficients(Runge_Kutta&);
            void Get_RK_VelocityContributions(double*, double*, double*);
            void Get_RK_ScalarContributions(double*, double*, double*);
            void Get_RK3_Integration(Mesher, Parallel);

            // Run Solver
            void RunSolver(Memory, ReadData, Parallel, Mesher, PostProcessing);
            
            // Memory Allocation for species
            void Allocate_StructSpecies(Memory, int);
            
            // Boundary Conditions for species
            void Get_InitialConditionsSpecies();
            void Get_InitialBoundaryConditionsSpecies();
            void Get_UpdateBoundaryConditionsSpecies(Mesher);
            void Get_InitialHalosSpecies();
            void Get_UpdateHalosSpecies();

            // Diffusion for species
            void Get_DiffusionCoefficients(int);
            void Get_DiffusionSpecies(Mesher, int);
            void Get_MassConservationSpecies(int);

            // Convection for species
            inline double ConvectiveSchemeSpecies(double, double, double, double, double, double, double, double, double, double, string);
            void Get_ConvectionSpecies(Mesher, int);

            void Get_Concentration();

};