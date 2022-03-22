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

    	int NX;
    	int NY;
    	int NZ;

		int Halo;
		int HP;
		
		double Tleft;
		double Tright;
		double Rayleigh;

		// Parameters for parallel computing
		int Rango;
		int Procesos;

		int *Ix;
		int *Fx;

		double *LocalNusselt;
		double *U_Averaged;
		double *V_Averaged;

		//Metodos de la clase
		void VTK_GlobalScalar3D(string, string, string, Mesher, double*);
		void VTK_GlobalVectorial3D(string, string, string, Mesher, double*, double*, double*);

		void Get_VelocityResults(Mesher, double*, double*, double);
		void Get_NusseltResults(Mesher, double*);

		void Delete_PostProcessingMemory();
};
