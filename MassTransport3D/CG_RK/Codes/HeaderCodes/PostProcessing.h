#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <mpi.h>

using namespace std;

class PostProcessing{	

	private:

	public:
		//Constructor de la clase
		PostProcessing(Memory, ReadData, Mesher, Parallel);
		
		//Datos de la clase
		int Problema;

		double Xdominio;
		double Ydominio;
		double Zdominio;

    	int NX;
    	int NY;
    	int NZ;

		int Halo;
		int HP;
		
		double Tleft;
		double Tright;
		double Rayleigh;
		double Reynolds;
		double Prandtl;

		// Parameters for parallel computing
		int Rango;
		int Procesos;

		int *Ix;
		int *Fx;

		double *LocalNusselt;
		double *U_Averaged;
		double *V_Averaged;

		double *Mass_Flow;
		double *Bulk_Temperature;
		double *Bulk_Fraction;
		double *T_Gradient;
		double *Y_Gradient;

		double *Channel_Nusselt;
		double *Channel_Sherwood;

		//Metodos de la clase
		void VTK_GlobalScalar3D(string, string, string, Mesher, double*);
		void VTK_GlobalVectorial3D(string, string, string, Mesher, double*, double*, double*);

		void Get_VelocityResults(Mesher, double*, double*, double);
		void Get_NusseltResults(Mesher, double*);

		void Get_NusseltNumber(Mesher, double*, double*, double*, double, double, double);
		void Get_SherwoodNumber(Mesher, double*, double*, double*, double);

		void Delete_PostProcessingMemory();
};
