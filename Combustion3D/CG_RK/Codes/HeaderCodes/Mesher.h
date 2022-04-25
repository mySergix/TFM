//------------------------------------------------------------------------------------------------//
//                             HEADER FILE FOR MESH CREATION CLASS                                //
//------------------------------------------------------------------------------------------------//

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

using namespace std;

class Mesher{
	private:

	public:
		//Constructor de la clase
		Mesher(Memory, ReadData, Parallel);

		// Datos Numéricos del problema
		string Problema;

		int NX;
		int NY; 
		int NZ;

		int NX_1, NX_2, NX_3;
		int NY_1, NY_2, NY_3;

		int OptionX_1, OptionX_2, OptionX_3;
		int OptionY_1, OptionY_2, OptionY_3;
		int OptionZ;

		double SFX_1, SFX_2, SFX_3;
		double SFY_1, SFY_2, SFY_3;
		double SFZ;

		// Datos Geométricos del problema
		double Width_Inlet;
		double Width_Slit;
		double Burner_Wall;
		double Symmetry_Burner;
		double Height_Slit;
		double Height_Burner;
		double Width_Burner;

		double X_1, X_2, X_3;
		double Y_1, Y_2, Y_3;

		double Xdomain, Ydomain, Zdomain;
		
		// Parámetros de computación paralela
		int Rango;
		int Procesos;
		int* Ix;
        int* Fx;

		int HP;
		int Halo;

		// Vector of columns
		double **NY_ColumnMP;
		double **NY_ColumnMU;
		double **NY_ColumnMV;
		double **NY_ColumnMW;

		// Meshes of the simulation
		double *MP; // Collocated mesh
		double *MU; // Staggered U mesh
		double *MV; // Staggered V mesh
		double *MW; // Staggered W mesh

		// Meshes nodal distances 
		double *DeltasMP; // Collocated mesh distances
		double *DeltasMU; // Staggered U mesh distances
		double *DeltasMV; // Staggered V mesh distances
		double *DeltasMW; // Staggered W mesh distances

		// Meshes surfaces of the CV
		double *SupMP; // Collocated mesh surfaces
		double *SupMU; // Staggered U mesh distances
		double *SupMV; // Staggered V mesh distances
		double *SupMW; // Staggered W mesh distances

		// Volumes of the CV of the meshes
		double *VolMP; // Collocated mesh CV Volumes
		double *VolMU; // Staggered U mesh CV Volumes
		double *VolMV; // Staggered V mesh CV Volumes
		double *VolMW; // Staggered W mesh CV Volumes


		// Global Mesh
		double *GlobalMeshP;
		double *GlobalMeshU;
		double *GlobalMeshV;
		double *GlobalMeshW;

		double *GlobalDeltasMV;
		
		void Allocate_MesherMemory(Memory); //Alojamiento de memoria para cada matriz
		void Delete_MesherMemory(); // Borrado de toda la memoria

		void Get_ColumnsNY();

		void Get_LocalMeshes(); //Creación de todas las mallas
		void Get_GlobalMesh(); // Creation of global collocated meesh

		void Get_Deltas(); //Cálculo de las distancias entre nodos en cada una de las matrices
		void Get_GlobalDeltas();

		void Get_Surfaces(); //Cálculo de las superficies de cada uno de los volúmenes de control
		void Get_Volumes(); //Cálculo de los volúmenes de control de cada volúmen
		
		//void MallaVTK3D(string, string, string, double*, int, int, int); //Pasar todos los resultados escalares en 2D a un archivo VTK
		void ExecuteMesher(Memory); //Ejecutar todos los procesos del mallador

};
		
