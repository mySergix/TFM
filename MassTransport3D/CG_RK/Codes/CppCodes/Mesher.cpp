//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR MESH CREATION CLASS                                   //
//------------------------------------------------------------------------------------------------//

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"

using namespace std;

//Constructor del mallador
Mesher::Mesher(Memory M1, ReadData R1, Parallel P1){
		
	//Datos Numéricos del problema
	Problema = R1.ProblemNumericalData[0];

	NX = R1.ProblemNumericalData[2];
	NY = R1.ProblemNumericalData[3];
	NZ = R1.ProblemNumericalData[4];
	Halo = 2;
	HP = 2;

	OptionX = R1.ProblemNumericalData[5];
	OptionY = R1.ProblemNumericalData[6];
	OptionZ = R1.ProblemNumericalData[7];

	SFX = R1.ProblemData[0];
	SFY = R1.ProblemData[1];
	SFZ = R1.ProblemData[2];

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
	
}

#include "Matrix_Index.cpp"
#include "Mesher_Memory.cpp"
#include "Mesher_NodalCoordinates.cpp"

//Cálculo de las distancias entre nodos en cada una de las matrices
void Mesher::Get_Deltas(){
int i, j, k;

	// Collocated mesh distances
	for(i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				DeltasMP[LP(i,j,k,0)] = MU[LU(i+1,j,k,0)] - MU[LU(i,j,k,0)]; // Delta X
				DeltasMP[LP(i,j,k,1)] = MV[LV(i,j+1,k,1)] - MV[LV(i,j,k,1)]; // Delta Y
				DeltasMP[LP(i,j,k,2)] = MW[LW(i,j,k+1,2)] - MW[LW(i,j,k,2)]; // Delta Z
			}
		}
	}

	// Staggered U mesh distances
	if (Rango != 0 && Rango != Procesos - 1){
		for(i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
			for(j = 0; j < NY; j++){
				for(k = 0; k < NZ; k++){
					DeltasMU[LU(i,j,k,0)] = MP[LP(i,j,k,0)] - MP[LP(i-1,j,k,0)]; // Delta X
					DeltasMU[LU(i,j,k,1)] = MV[LV(i,j+1,k,1)] - MV[LV(i,j,k,1)]; // Delta Y
					DeltasMU[LU(i,j,k,2)] = MW[LW(i,j,k+1,2)] - MW[LW(i,j,k,2)]; // Delta Z
				}
			}
		}
	}
	else if (Rango == 0){ // Left
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){

				DeltasMU[LU(0,j,k,0)] = MP[LP(0,j,k,0)] - MU[LU(0,j,k,0)]; // Delta X Parte Izquierda
				DeltasMU[LU(0,j,k,1)] = MV[LV(0,j+1,k,1)] - MV[LV(0,j,k,1)]; // Delta Y Parte Izquierda
				DeltasMU[LU(0,j,k,2)] = MW[LW(0,j,k+1,2)] - MW[LW(0,j,k,2)]; // Delta Z Parte Izquierda

				for(i = Ix[Rango] + 1; i < Fx[Rango] + 1; i++){
					DeltasMU[LU(i,j,k,0)] = MP[LP(i,j,k,0)] - MP[LP(i-1,j,k,0)]; // Delta X
					DeltasMU[LU(i,j,k,1)] = MV[LV(i,j+1,k,1)] - MV[LV(i,j,k,1)]; // Delta Y
					DeltasMU[LU(i,j,k,2)] = MW[LW(i,j,k+1,2)] - MW[LW(i,j,k,2)]; // Delta Z
				}
			}
		}
	}
	else if(Rango == Procesos - 1){ // Right
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){

				DeltasMU[LU(NX,j,k,0)] = MU[LU(NX,j,k,0)] - MP[LP(NX-1,j,k,0)]; //Delta X Parte Derecha
				DeltasMU[LU(NX,j,k,1)] = MV[LV(NX-1,j+1,k,1)] - MV[LV(NX-1,j,k,1)]; //Delta Y Parte Derecha
				DeltasMU[LU(NX,j,k,2)] = MW[LW(NX-1,j,k+1,2)] - MW[LW(NX-1,j,k,2)]; //Delta Z Parte Derecha

				for(i = Ix[Rango] - 1; i < Fx[Rango]; i++){
					DeltasMU[LU(i,j,k,0)] = MP[LP(i,j,k,0)] - MP[LP(i-1,j,k,0)]; // Delta X
					DeltasMU[LU(i,j,k,1)] = MV[LV(i,j+1,k,1)] - MV[LV(i,j,k,1)]; // Delta Y
					DeltasMU[LU(i,j,k,2)] = MW[LW(i,j,k+1,2)] - MW[LW(i,j,k,2)]; // Delta Z
				}
			}
		}
	}

	// Staggered V mesh distances
	for(i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		for(k = 0; k < NZ; k++){

			DeltasMV[LV(i,0,k,0)] = MU[LU(i+1,0,k,0)] - MU[LU(i,0,k,0)]; // Delta X Parte Abajo
			DeltasMV[LV(i,NY,k,0)] = MU[LU(i+1,NY-1,k,0)] - MU[LU(i,NY-1,k,0)]; // Delta X Parte Arriba

			DeltasMV[LV(i,0,k,1)] = MP[LP(i,0,k,1)] - MV[LV(i,0,k,1)]; // Delta Y Parte Abajo
			DeltasMV[LV(i,NY,k,1)] = MV[LV(i,NY,k,1)] - MP[LP(i,NY-1,k,1)]; // Delta Y Parte Arriba

			DeltasMV[LV(i,0,k,2)] = MW[LW(i,0,k+1,2)] - MW[LW(i,0,k,2)]; // Delta Z Parte Abajo
			DeltasMV[LV(i,NY,k,2)] = MW[LW(i,NY-1,k+1,2)] - MW[LW(i,NY-1,k,2)]; // Delta Z Parte Arriba

			for(j = 1; j < NY; j++){
				DeltasMV[LV(i,j,k,0)] = MU[LU(i+1,j,k,0)] - MU[LU(i,j,k,0)]; //Delta X
				DeltasMV[LV(i,j,k,1)] = MP[LP(i,j,k,1)] - MP[LP(i,j-1,k,1)]; //Delta Y
				DeltasMV[LV(i,j,k,2)] = MW[LW(i,j,k+1,2)] - MW[LW(i,j,k,2)]; //Delta Z
			}

		}
	}

	// Staggered W mesh distances
	for(i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){

			DeltasMW[LW(i,j,0,0)] = MU[LU(i+1,j,0,0)] - MU[LU(i,j,0,0)]; // Delta X Parte Here
			DeltasMW[LW(i,j,NZ,0)] = MU[LU(i+1,j,NZ-1,0)] - MU[LU(i,j,NZ-1,0)]; // Delta X Parte There

			DeltasMW[LW(i,j,0,1)] = MV[LV(i,j+1,0,1)] - MV[LV(i,j,0,1)]; // Delta Y Parte Here
			DeltasMW[LW(i,j,NZ,1)] = MV[LV(i,j+1,NZ-1,1)] - MV[LV(i,j,NZ-1,1)]; // Delta Y Parte There

			DeltasMW[LW(i,j,0,2)] = MP[LP(i,j,0,2)] - MW[LW(i,j,0,2)]; // Delta Z Parte Here
			DeltasMW[LW(i,j,NZ,2)] = MW[LW(i,j,NZ,2)] - MP[LP(i,j,NZ-1,2)]; // Delta Z Parte There

			for(k = 1; k < NZ; k++){
				DeltasMW[LW(i,j,k,0)] = MU[LU(i+1,j,k,0)] - MU[LU(i,j,k,0)]; //Delta X
				DeltasMW[LW(i,j,k,1)] = MV[LV(i,j+1,k,1)] - MV[LV(i,j,k,1)]; //Delta Y
				DeltasMW[LW(i,j,k,2)] = MP[LP(i,j,k,2)] - MP[LP(i,j,k-1,2)]; //Delta Z
			}
		}
	}

}

//Cálculo de las superficies de cada uno de los volúmenes de control
void Mesher::Get_Surfaces(){
int i, j, k;

	// Collocated mesh surfaces
	for(i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				SupMP[LP(i,j,k,0)] = DeltasMU[LU(i,j,k,1)]*DeltasMU[LU(i,j,k,2)]; // Surface X
				SupMP[LP(i,j,k,1)] = DeltasMV[LV(i,j,k,0)]*DeltasMV[LV(i,j,k,2)]; // Surface Y
				SupMP[LP(i,j,k,2)] = DeltasMW[LW(i,j,k,0)]*DeltasMW[LW(i,j,k,1)]; // Surface Z
			}
		}
	}

	// Staggered U mesh surfaces
	for(i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				SupMU[LU(i,j,k,0)] = DeltasMU[LU(i,j,k,1)]*DeltasMU[LU(i,j,k,2)]; // Surface X
				SupMU[LU(i,j,k,1)] = DeltasMU[LU(i,j,k,0)]*DeltasMU[LU(i,j,k,2)]; // Surface Y 
				SupMU[LU(i,j,k,2)] = DeltasMU[LU(i,j,k,0)]*DeltasMU[LU(i,j,k,1)]; // Surface Z
			}
		}
	}	

	// Staggered V mesh surfaces
	for(i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY + 1; j++){
				SupMV[LV(i,j,k,0)] = DeltasMV[LV(i,j,k,1)]*DeltasMV[LV(i,j,k,2)]; // Surface X
				SupMV[LV(i,j,k,1)] = DeltasMV[LV(i,j,k,0)]*DeltasMV[LV(i,j,k,2)]; // Surface Y 
				SupMV[LV(i,j,k,2)] = DeltasMV[LV(i,j,k,0)]*DeltasMV[LV(i,j,k,1)]; // Surface Z
			}
		}
	}

	// Staggered W mesh surfaces
	for(i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){		
			for(k = 0; k < NZ + 1; k++){
				SupMW[LW(i,j,k,0)] = DeltasMW[LW(i,j,k,1)]*DeltasMW[LW(i,j,k,2)]; // Surface X
				SupMW[LW(i,j,k,1)] = DeltasMW[LW(i,j,k,0)]*DeltasMW[LW(i,j,k,2)]; // Surface Y 
				SupMW[LW(i,j,k,2)] = DeltasMW[LW(i,j,k,0)]*DeltasMW[LW(i,j,k,1)]; // Surface Z
			}
		}
	}

}

//Cálculo de los volúmenes de control de cada volúmen
void Mesher::Get_Volumes(){
int i, j, k;

	// Collocated mesh CV Volumes
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				VolMP[LP(i,j,k,0)] = DeltasMP[LP(i,j,k,0)]*DeltasMP[LP(i,j,k,1)]*DeltasMP[LP(i,j,k,2)];
			}
		}
	}

	// Staggered U mesh CV Volumes
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				VolMU[LU(i,j,k,0)] = DeltasMU[LU(i,j,k,0)]*DeltasMU[LU(i,j,k,1)]*DeltasMU[LU(i,j,k,2)];
			}
		}
	}

	// Staggered V mesh CV Volumes
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY + 1; j++){
			for(k = 0; k < NZ; k++){
				VolMV[LV(i,j,k,0)] = DeltasMV[LV(i,j,k,0)]*DeltasMV[LV(i,j,k,1)]*DeltasMV[LV(i,j,k,2)];
			}
		}
	}

	// Staggered W mesh CV Volumes
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ + 1; k++){
				VolMW[LW(i,j,k,0)] = DeltasMW[LW(i,j,k,0)]*DeltasMW[LW(i,j,k,1)]*DeltasMW[LW(i,j,k,2)];
			}
		}
	}

}

// Function to calulate the global volume distances
void Mesher::Get_GlobalDeltas(){
int i, j, k;

	// Collocated mesh distances
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				GlobalDeltasMP[GP(i,j,k,0)] = GlobalMeshU[GU(i+1,j,k,0)] - GlobalMeshU[GU(i,j,k,0)]; // Delta X
				GlobalDeltasMP[GP(i,j,k,1)] = GlobalMeshV[GV(i,j+1,k,1)] - GlobalMeshV[GV(i,j,k,1)]; // Delta Y
				GlobalDeltasMP[GP(i,j,k,2)] = GlobalMeshW[GW(i,j,k+1,2)] - GlobalMeshW[GW(i,j,k,2)]; // Delta Z
			}
		}
	}

	// Staggered U mesh distances
	for(i = 1; i < NX; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				GlobalDeltasMU[GU(i,j,k,0)] = GlobalMeshP[GP(i,j,k,0)] - GlobalMeshP[GP(i-1,j,k,0)]; // Delta X
				GlobalDeltasMU[GU(i,j,k,1)] = GlobalMeshV[GV(i,j+1,k,1)] - GlobalMeshV[GV(i,j,k,1)]; // Delta Y
				GlobalDeltasMU[GU(i,j,k,2)] = GlobalMeshW[GW(i,j,k+1,2)] - GlobalMeshW[GW(i,j,k,2)]; // Delta Z
			}
		}
	}
	
	for(j = 0; j < NY; j++){
		for(k = 0; k < NZ; k++){
			GlobalDeltasMU[GU(0,j,k,0)] = GlobalMeshP[GP(0,j,k,0)] - GlobalMeshU[GU(0,j,k,0)]; // Delta X Parte Izquierda
			GlobalDeltasMU[GU(0,j,k,1)] = GlobalMeshV[GV(0,j+1,k,1)] - GlobalMeshV[GV(0,j,k,1)]; // Delta Y Parte Izquierda
			GlobalDeltasMU[GU(0,j,k,2)] = GlobalMeshW[GW(0,j,k+1,2)] - GlobalMeshW[GW(0,j,k,2)]; // Delta Z Parte Izquierda	
		}
	}

	for(j = 0; j < NY; j++){
		for(k = 0; k < NZ; k++){
			GlobalDeltasMU[GU(NX,j,k,0)] = GlobalMeshU[GU(NX,j,k,0)] - GlobalMeshP[GP(NX-1,j,k,0)]; //Delta X Parte Derecha
			GlobalDeltasMU[GU(NX,j,k,1)] = GlobalMeshV[GV(NX-1,j+1,k,1)] - GlobalMeshV[GV(NX-1,j,k,1)]; //Delta Y Parte Derecha
			GlobalDeltasMU[GU(NX,j,k,2)] = GlobalMeshW[GW(NX-1,j,k+1,2)] - GlobalMeshW[GW(NX-1,j,k,2)]; //Delta Z Parte Derecha			
		}
	}

	// Staggered V mesh distances
	for(i = 0; i < NX; i++){
		for(k = 0; k < NZ; k++){

			GlobalDeltasMV[GV(i,0,k,0)] = GlobalMeshU[GU(i+1,0,k,0)] - GlobalMeshU[GU(i,0,k,0)]; // Delta X Parte Abajo
			GlobalDeltasMV[GV(i,NY,k,0)] = GlobalMeshU[GU(i+1,NY-1,k,0)] - GlobalMeshU[GU(i,NY-1,k,0)]; // Delta X Parte Arriba

			GlobalDeltasMV[GV(i,0,k,1)] = GlobalMeshP[GP(i,0,k,1)] - GlobalMeshV[GV(i,0,k,1)]; // Delta Y Parte Abajo
			GlobalDeltasMV[GV(i,NY,k,1)] = GlobalMeshV[GV(i,NY,k,1)] - GlobalMeshP[GP(i,NY-1,k,1)]; // Delta Y Parte Arriba

			GlobalDeltasMV[GV(i,0,k,2)] = GlobalMeshW[GW(i,0,k+1,2)] - GlobalMeshW[GW(i,0,k,2)]; // Delta Z Parte Abajo
			GlobalDeltasMV[GV(i,NY,k,2)] = GlobalMeshW[GW(i,NY-1,k+1,2)] - GlobalMeshW[GW(i,NY-1,k,2)]; // Delta Z Parte Arriba

			for(j = 1; j < NY; j++){
				GlobalDeltasMV[GV(i,j,k,0)] = GlobalMeshU[GU(i+1,j,k,0)] - GlobalMeshU[GU(i,j,k,0)]; //Delta X
				GlobalDeltasMV[GV(i,j,k,1)] = GlobalMeshP[GP(i,j,k,1)] - GlobalMeshP[GP(i,j-1,k,1)]; //Delta Y
				GlobalDeltasMV[GV(i,j,k,2)] = GlobalMeshW[GW(i,j,k+1,2)] - GlobalMeshW[GW(i,j,k,2)]; //Delta Z
			}

		}
	}

	// Staggered W mesh distances
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){

			GlobalDeltasMW[GW(i,j,0,0)] = GlobalMeshU[GU(i+1,j,0,0)] - GlobalMeshU[GU(i,j,0,0)]; // Delta X Parte Here
			GlobalDeltasMW[GW(i,j,NZ,0)] = GlobalMeshU[GU(i+1,j,NZ-1,0)] - GlobalMeshU[GU(i,j,NZ-1,0)]; // Delta X Parte There

			GlobalDeltasMW[GW(i,j,0,1)] = GlobalMeshV[GV(i,j+1,0,1)] - GlobalMeshV[GV(i,j,0,1)]; // Delta Y Parte Here
			GlobalDeltasMW[GW(i,j,NZ,1)] = GlobalMeshV[GV(i,j+1,NZ-1,1)] - GlobalMeshV[GV(i,j,NZ-1,1)]; // Delta Y Parte There

			GlobalDeltasMW[GW(i,j,0,2)] = GlobalMeshP[GP(i,j,0,2)] - GlobalMeshW[GW(i,j,0,2)]; // Delta Z Parte Here
			GlobalDeltasMW[GW(i,j,NZ,2)] = GlobalMeshW[GW(i,j,NZ,2)] - GlobalMeshP[GP(i,j,NZ-1,2)]; // Delta Z Parte There

			for(k = 1; k < NZ; k++){
				GlobalDeltasMW[GW(i,j,k,0)] = GlobalMeshU[GU(i+1,j,k,0)] - GlobalMeshU[GU(i,j,k,0)]; //Delta X
				GlobalDeltasMW[GW(i,j,k,1)] = GlobalMeshV[GV(i,j+1,k,1)] - GlobalMeshV[GV(i,j,k,1)]; //Delta Y
				GlobalDeltasMW[GW(i,j,k,2)] = GlobalMeshP[GP(i,j,k,2)] - GlobalMeshP[GP(i,j,k-1,2)]; //Delta Z
			}
		}
	}

}

// Function to calculate the global mesh surfaces
void Mesher::Get_GlobalSurfaces(){
int i, j, k;

	// Collocated mesh surfaces
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				GlobalSupMP[GP(i,j,k,0)] = GlobalDeltasMU[GU(i,j,k,1)]*GlobalDeltasMU[GU(i,j,k,2)]; // Surface X
				GlobalSupMP[GP(i,j,k,1)] = GlobalDeltasMV[GV(i,j,k,0)]*GlobalDeltasMV[GV(i,j,k,2)]; // Surface Y
				GlobalSupMP[GP(i,j,k,2)] = GlobalDeltasMW[GW(i,j,k,0)]*GlobalDeltasMW[GW(i,j,k,1)]; // Surface Z
			}
		}
	}

}
//Ejecutar todos los procesos del mallador
void Mesher::ExecuteMesher(Memory M1){
	
	Allocate_MesherMemory(M1); // Memory Allocation
	Get_LocalMeshes(); // Creación de todas las mallas
	Get_Deltas(); // Cálculo de las distancias entre nodos en cada una de las matrices
	Get_Surfaces(); // Cálculo de las superficies de cada uno de los volúmenes de control
	Get_Volumes(); // Cálculo de los volúmenes de control de cada volúmen
	
	if(Rango == 0){	
		Get_GlobalMesh();
		Get_GlobalDeltas();
		Get_GlobalSurfaces();
	//	MallaVTK3D("ParaviewResults/MeshResults/", "MallaP", "MallaMP", MP, NX, NY, NZ);
		cout<<"Mesh created."<<endl;
	}

}

