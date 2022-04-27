//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR PARALLEL PROGRAMMING CLASS                            //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"

using namespace std;

Parallel::Parallel(ReadData R1){

}

// Additional files of the class
#include "Matrix_Index.cpp"
#include "Parallel_Communication.cpp"

// Function to get the rank of each core
void Parallel::Rango_Procesos(){
	int a;
	a = MPI_Comm_rank(MPI_COMM_WORLD, &Rango);
}

// Function to get the total number of cores
void Parallel::Total_Procesos(){
	int a;
	a = MPI_Comm_size(MPI_COMM_WORLD, &Procesos);
}

// Function to allocate memory for all cores Ix and Fx
void Parallel::AllocateMemory(Memory M1){
    Ix = M1.AllocateInt(Procesos, 1, 1, 1);
    Fx = M1.AllocateInt(Procesos, 1, 1, 1);
}

// Function to perform the worksplit between all cores
void Parallel::WorkSplit(int NX, int *Ix, int *Fx){
MPI_Status ST;
int i;
int Intervalo, Residuo;
int p;

	Intervalo = NX/Procesos;
	Residuo = NX%Procesos;

	if(Rango != Procesos-1){
		Ix[Rango] = Rango*Intervalo;
		Fx[Rango] = (Rango+1)*Intervalo;
	}
	else{
		Ix[Rango] = Rango*Intervalo;
		Fx[Rango] = (Rango+1)*Intervalo + Residuo;
	}

    // Send Ix values to every processor
    for (i = 0; i < Procesos; i++){
        if (i != Rango){
            MPI_Send(&Ix[Rango], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i < Procesos; i++){
        if(i != Rango){
            MPI_Recv(&Ix[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ST);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Send Fx values to every processor
    for (i = 0; i < Procesos; i++){
        if (i != Rango){
            MPI_Send(&Fx[Rango], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i < Procesos; i++){
        if(i != Rango){
            MPI_Recv(&Fx[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &ST);
        }
    }

}

// Function to execute the parallel class
void Parallel::RunParallel(Memory M1){

    Rango_Procesos();
    Total_Procesos();
    AllocateMemory(M1);
    
}

// Function to get the information of the mesh splits
void Parallel::Get_MesherInformation(Memory M1, int *Mesh_Ix, int *Mesh_Fx, int Mesh_NX, int Mesh_NY, int Mesh_NZ, int Mesh_HP, int Mesh_Halo, int **NY_MP, int **NY_MU, int **NY_MV, int **NY_MW){
int i;

    for (i = 0; i < Procesos; i++){
        Ix[i] = Mesh_Ix[i];
        Fx[i] = Mesh_Fx[i];
    }

    NX = Mesh_NX;
    NY = Mesh_NY;
    NZ = Mesh_NZ;

    HP = Mesh_HP;
    Halo = Mesh_Halo;

    // Vector of columns Memory Allocation
	NY_ColumnMP = M1.AllocateInt_Matrix2D(Fx[Rango] - Ix[Rango] + 2*HP, 2);
	NY_ColumnMU = M1.AllocateInt_Matrix2D(Fx[Rango] - Ix[Rango] + 2*Halo + 1, 2);
	NY_ColumnMV = M1.AllocateInt_Matrix2D(Fx[Rango] - Ix[Rango] + 2*Halo, 2);
	NY_ColumnMW = M1.AllocateInt_Matrix2D(Fx[Rango] - Ix[Rango] + 2*Halo, 2);

    // Staggered Mesh U Columns
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo + 1; i++){
        NY_ColumnMU[i + Halo - Ix[Rango]][0] = NY_MU[i + Halo - Ix[Rango]][0]; // Bottom
        NY_ColumnMU[i + Halo - Ix[Rango]][1] = NY_MU[i + Halo - Ix[Rango]][1]; // Top
    }

    // Meshes P, V and W Columns
    for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
        NY_ColumnMP[i + Halo - Ix[Rango]][0] = NY_MP[i + Halo - Ix[Rango]][0]; // Bottom
        NY_ColumnMP[i + Halo - Ix[Rango]][1] = NY_MP[i + Halo - Ix[Rango]][1]; // Top

        NY_ColumnMV[i + Halo - Ix[Rango]][0] = NY_MV[i + Halo - Ix[Rango]][0]; // Bottom
        NY_ColumnMV[i + Halo - Ix[Rango]][1] = NY_MV[i + Halo - Ix[Rango]][1]; // Top

        NY_ColumnMW[i + Halo - Ix[Rango]][0] = NY_MW[i + Halo - Ix[Rango]][0]; // Bottom
        NY_ColumnMW[i + Halo - Ix[Rango]][1] = NY_MW[i + Halo - Ix[Rango]][1]; // Top
    }

}

// Function to delete all the Parallel class memory
void Parallel::Delete_ParallelMemory(){

    delete[] Ix;
    delete[] Fx;
    
    delete[] NY_ColumnMP;
    delete[] NY_ColumnMU;
    delete[] NY_ColumnMV;
    delete[] NY_ColumnMW;
    
}