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
		string Problema;

    	int NX;
    	int NY;
    	int NZ;

		int Halo;
		int HP;
		
		// Parameters for parallel computing
		int Rango;
		int Procesos;

		int *Ix;
		int *Fx;

		//Metodos de la clase
		void Get_GlobalVectorialHalos(double*, double*, double*);

		void VTK_GlobalScalar3D(string, string, string, Mesher, double*);
		void VTK_GlobalVectorial3D(string, string, string, Mesher, double*, double*, double*);

		void Delete_PostProcessingMemory();
};
