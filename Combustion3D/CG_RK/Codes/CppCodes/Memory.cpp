//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR MEMORY ALLOCATION CLASS                               //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"

using namespace std;

Memory::Memory(){
	
};

//Memoria dinámica (double)
double *Memory::AllocateDouble(int NX, int NY, int NZ, int Dim){
double *M1;

	M1 = new double [NX*NY*NZ*Dim];				
	return M1;

}

//Memoria dinámica (int)
int *Memory::AllocateInt(int NX, int NY, int NZ, int Dim){ 
int *M1;

	M1 = new int [NX*NY*NZ*Dim];				
	return M1;

}

// Memoria dinámica matriz (int)
int **Memory::AllocateInt_Matrix2D(int NX, int NY){ 
int i;
int **M1;

	M1 = new int *[NX];		

	for (i = 0; i < NX; i++){
		M1[i] = new int [NY];	
	}	

	return M1;

}



