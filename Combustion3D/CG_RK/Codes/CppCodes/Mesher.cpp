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
	Problema = "Premixed";

	NX_1 = 20;
	NX_2 = 10;
	NX_3 = 30;

	NY_1 = 30;
	NY_2 = 10;
	NY_3 = 40;
	NY_4 = 50;

	// Meshing options:
	// 1 -> Regular

		// OptionX:
			// 2 -> Right Sided - Hyperbolic Tangent (More Density on the Right)
			// 3 -> Left Sided - Hyperbolic Tangent (More Density on the Left)
			// 4 -> Centered Sided - Hyperbolic Tangent (More Density on both Sides)

		// OptionY:
			// 2 -> Top Sided - Hyperbolic Tangent (More Density on the Top)
			// 3 -> Down Sided - Hyperbolic Tangent (More Density on the Down)
			// 4 -> Centered Sided - Hyperbolic Tangent (More Density on both Sides)

	OptionX_1 = 1;
	OptionX_2 = 1;
	OptionX_3 = 1;

	OptionY_1 = 1;
	OptionY_2 = 1;
	OptionY_3 = 1;
	OptionY_4 = 1;

	OptionZ = 1;

	SFX_1 = 2.0;
	SFX_2 = 2.0;
	SFX_3 = 2.0;

	SFY_1 = 2.0;
	SFY_2 = 2.0;
	SFY_3 = 2.0;
	SFY_4 = 2.0;

	SFZ = 1.0;

	NX = NX_1 + NX_2 + NX_3;
	if (Problema == "Premixed"){ NY = NY_1 + NY_2 + NY_3 + NY_4; }
	else if (Problema == "NonPremixed"){ NY = NY_2 + NY_3 + NY_4; }
	NZ = 10;

	// Datos Geométricos del problema (m)
	Width_Inlet = 0.002; 
	Width_Slit = 0.001;
	Burner_Wall = 0.003;
	Symmetry_Burner = 0.013;
	Height_Slit = 0.001;
	Height_Burner =  0.010;
	Width_Burner = 0.006;

	FlameFront = 0.004;

	X_1 = Width_Inlet;
	X_2 = Width_Slit;
	X_3 = Burner_Wall;

	Y_1 = Symmetry_Burner - Height_Burner;
	Y_2 = Height_Slit;
	Y_3 = FlameFront;
	Y_4 = Height_Burner - Height_Slit - FlameFront;

	Xdomain = X_1 + X_2 + X_3;
	if (Problema == "Premixed"){ Ydomain = Y_1 + Y_2 + Y_3 + Y_4; }
	else if (Problema == "NonPremixed"){ Ydomain = Y_2 + Y_3 + Y_4; }
	Zdomain = 0.002;

	//Datos necesarios para computación paralela
	Rango = P1.Rango;
	Procesos = P1.Procesos;
	Ix = M1.AllocateInt(Procesos, 1, 1, 1);
    Fx = M1.AllocateInt(Procesos, 1, 1, 1);

	P1.WorkSplit(NX, Ix, Fx);

	Halo = 2;
	HP = 2;

}

#include "Matrix_Index.cpp"
#include "Mesher_Memory.cpp"
#include "Mesher_NodalCoordinates.cpp"
#include "Mesher_Paraview.cpp"

// Seteo del numero de nodos NY por cada columna en direccion X (Local)
void Mesher::Get_LocalColumnsNY(){
int i, j;

	if (Problema == "Premixed"){

		// Nodes U
		for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo + 1; i++){

			if (i <= NX_1){
				NY_ColumnMU[i + Halo - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMU[i + Halo - Ix[Rango]][1] = NY; // Top
			}
			else if (i > NX_1 && i < NX_1 + NX_2){
				NY_ColumnMU[i + Halo - Ix[Rango]][0] = NY_1 + NY_2; // Bottom
				NY_ColumnMU[i + Halo - Ix[Rango]][1] = NY; // Top
			}
			else if (i >= NX_1 + NX_2){
				NY_ColumnMU[i + Halo - Ix[Rango]][0] = NY_1; // Bottom
				NY_ColumnMU[i + Halo - Ix[Rango]][1] = NY; // Top
			}

		}

		// Nodes P, V and W
		for (i = Ix[Rango] - HP; i < Fx[Rango] + HP; i++){

			if (i < NX_1){
				NY_ColumnMP[i + HP - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMP[i + HP - Ix[Rango]][1] = NY; // Top

				NY_ColumnMV[i + Halo - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMV[i + Halo - Ix[Rango]][1] = NY + 1; // Top

				NY_ColumnMW[i + Halo - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMW[i + Halo - Ix[Rango]][1] = NY; // Top
			}
			else if (i >= NX_1 && i < NX_1 + NX_2){
				NY_ColumnMP[i + HP - Ix[Rango]][0] = NY_1 + NY_2; // Bottom
				NY_ColumnMP[i + HP - Ix[Rango]][1] = NY; // Top

				NY_ColumnMV[i + Halo - Ix[Rango]][0] = NY_1 + NY_2; // Bottom
				NY_ColumnMV[i + Halo - Ix[Rango]][1] = NY + 1; // Top

				NY_ColumnMW[i + Halo - Ix[Rango]][0] = NY_1 + NY_2; // Bottom
				NY_ColumnMW[i + Halo - Ix[Rango]][1] = NY; // Top
			}
			else if (i >= NX_1 + NX_2){
				NY_ColumnMP[i + HP - Ix[Rango]][0] = NY_1; // Bottom
				NY_ColumnMP[i + HP - Ix[Rango]][1] = NY; // Top

				NY_ColumnMV[i + Halo - Ix[Rango]][0] = NY_1; // Bottom
				NY_ColumnMV[i + Halo - Ix[Rango]][1] = NY + 1; // Top

				NY_ColumnMW[i + Halo - Ix[Rango]][0] = NY_1; // Bottom
				NY_ColumnMW[i + Halo - Ix[Rango]][1] = NY; // Top
			}

		}

	}
	else if (Problema == "NonPremixed"){

		// Nodes U
		for (i = Ix[Rango] - HP; i < Fx[Rango] + HP; i++){

			if (i <= NX_1){
				NY_ColumnMU[i + Halo - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMU[i + Halo - Ix[Rango]][1] = NY; // Top
			}
			else if (i > NX_1 && i < NX_1 + NX_2){
				NY_ColumnMU[i + Halo - Ix[Rango]][0] = NY_2; // Bottom
				NY_ColumnMU[i + Halo - Ix[Rango]][1] = NY; // Top
			}
			else if (i >= NX_1 + NX_2){
				NY_ColumnMU[i + Halo - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMU[i + Halo - Ix[Rango]][1] = NY; // Top
			}

		}
		
		// Nodes P, V and W
		for (i = Ix[Rango] - HP; i < Fx[Rango] + HP; i++){

			if (i < NX_1){
				NY_ColumnMP[i + HP - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMP[i + HP - Ix[Rango]][1] = NY; // Top

				NY_ColumnMV[i + Halo - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMV[i + Halo - Ix[Rango]][1] = NY + 1; // Top

				NY_ColumnMW[i + Halo - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMW[i + Halo - Ix[Rango]][1] = NY; // Top
			}
			else if (i >= NX_1 && i < NX_1 + NX_2){
				NY_ColumnMP[i + HP - Ix[Rango]][0] = NY_2; // Bottom
				NY_ColumnMP[i + HP - Ix[Rango]][1] = NY; // Top

				NY_ColumnMV[i + Halo - Ix[Rango]][0] = NY_2; // Bottom
				NY_ColumnMV[i + Halo - Ix[Rango]][1] = NY + 1; // Top

				NY_ColumnMW[i + Halo - Ix[Rango]][0] = NY_2; // Bottom
				NY_ColumnMW[i + Halo - Ix[Rango]][1] = NY; // Top
			}
			else if (i >= NX_1 + NX_2){
				NY_ColumnMP[i + HP - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMP[i + HP - Ix[Rango]][1] = NY; // Top

				NY_ColumnMV[i + Halo - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMV[i + Halo - Ix[Rango]][1] = NY + 1; // Top

				NY_ColumnMW[i + Halo - Ix[Rango]][0] = 0; // Bottom
				NY_ColumnMW[i + Halo - Ix[Rango]][1] = NY; // Top
			}

		}

	}
	
}

// Seteo del numero de nodos NY por cada columna en direccion X (Global)
void Mesher::Get_GlobalColumnsNY(){
int i;

	if (Problema == "Premixed"){

		// Nodes U
		for (i = 0; i < NX; i++){

			if (i <= NX_1){
				GlobalNY_ColumnMU[i][0] = 0; // Bottom
				GlobalNY_ColumnMU[i][1] = NY; // Top
			}
			else if (i > NX_1 && i < NX_1 + NX_2){
				GlobalNY_ColumnMU[i][0] = NY_1 + NY_2; // Bottom
				GlobalNY_ColumnMU[i][1] = NY; // Top
			}
			else if (i >= NX_1 + NX_2){
				GlobalNY_ColumnMU[i][0] = NY_1; // Bottom
				GlobalNY_ColumnMU[i][1] = NY; // Top
			}

		}

		// Nodes P, V and W
		for (i = 0; i < NX; i++){

			if (i < NX_1){
				GlobalNY_ColumnMP[i][0] = 0; // Bottom
				GlobalNY_ColumnMP[i][1] = NY; // Top

				GlobalNY_ColumnMV[i][0] = 0; // Bottom
				GlobalNY_ColumnMV[i][1] = NY + 1; // Top

				GlobalNY_ColumnMW[i][0] = 0; // Bottom
				GlobalNY_ColumnMW[i][1] = NY; // Top
			}
			else if (i >= NX_1 && i < NX_1 + NX_2){
				GlobalNY_ColumnMP[i][0] = NY_1 + NY_2; // Bottom
				GlobalNY_ColumnMP[i][1] = NY; // Top

				GlobalNY_ColumnMV[i][0] = NY_1 + NY_2; // Bottom
				GlobalNY_ColumnMV[i][1] = NY + 1; // Top

				GlobalNY_ColumnMW[i][0] = NY_1 + NY_2; // Bottom
				GlobalNY_ColumnMW[i][1] = NY; // Top
			}
			else if (i >= NX_1 + NX_2){
				GlobalNY_ColumnMP[i][0] = NY_1; // Bottom
				GlobalNY_ColumnMP[i][1] = NY; // Top

				GlobalNY_ColumnMV[i][0] = NY_1; // Bottom
				GlobalNY_ColumnMV[i][1] = NY + 1; // Top

				GlobalNY_ColumnMW[i][0] = NY_1; // Bottom
				GlobalNY_ColumnMW[i][1] = NY; // Top
			}

		}

	}
	else if (Problema == "NonPremixed"){

		// Nodes U
		for (i = 0; i < NX; i++){

			if (i <= NX_1){
				GlobalNY_ColumnMU[i][0] = 0; // Bottom
				GlobalNY_ColumnMU[i][1] = NY; // Top
			}
			else if (i > NX_1 && i < NX_1 + NX_2){
				GlobalNY_ColumnMU[i][0] = NY_2; // Bottom
				GlobalNY_ColumnMU[i][1] = NY; // Top
			}
			else if (i >= NX_1 + NX_2){
				GlobalNY_ColumnMU[i][0] = 0; // Bottom
				GlobalNY_ColumnMU[i][1] = NY; // Top
			}

		}
		
		// Nodes P, V and W
		for (i = 0; i < NX; i++){

			if (i < NX_1){
				GlobalNY_ColumnMP[i][0] = 0; // Bottom
				GlobalNY_ColumnMP[i][1] = NY; // Top

				GlobalNY_ColumnMV[i][0] = 0; // Bottom
				GlobalNY_ColumnMV[i][1] = NY + 1; // Top

				GlobalNY_ColumnMW[i][0] = 0; // Bottom
				GlobalNY_ColumnMW[i][1] = NY; // Top
			}
			else if (i >= NX_1 && i < NX_1 + NX_2){
				GlobalNY_ColumnMP[i][0] = NY_2; // Bottom
				GlobalNY_ColumnMP[i][1] = NY; // Top

				GlobalNY_ColumnMV[i][0] = NY_2; // Bottom
				GlobalNY_ColumnMV[i][1] = NY + 1; // Top

				GlobalNY_ColumnMW[i][0] = NY_2; // Bottom
				GlobalNY_ColumnMW[i][1] = NY; // Top
			}
			else if (i >= NX_1 + NX_2){
				GlobalNY_ColumnMP[i][0] = 0; // Bottom
				GlobalNY_ColumnMP[i][1] = NY; // Top

				GlobalNY_ColumnMV[i][0] = 0; // Bottom
				GlobalNY_ColumnMV[i][1] = NY + 1; // Top

				GlobalNY_ColumnMW[i][0] = 0; // Bottom
				GlobalNY_ColumnMW[i][1] = NY; // Top
			}

		}

	}
	
}

// Function to calculate the total number of nodes in the global mesh
void Mesher::Get_TotalNodes(){
int i;
TotalNodesP = 0;

	// Collocated Pressure Nodes
	for (i = 0; i < NX; i++){
		TotalNodesP += (GlobalNY_ColumnMP[i][1] - GlobalNY_ColumnMP[i][0]) * NZ;
	}

}

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
				GlobalDeltasMP[GP(i,j,k,0)] = GlobalMeshU[LU(i+1,j,k,0)] - GlobalMeshU[GU(i,j,k,0)]; // Delta X
				GlobalDeltasMP[GP(i,j,k,1)] = GlobalMeshV[LV(i,j+1,k,1)] - GlobalMeshV[GV(i,j,k,1)]; // Delta Y
				GlobalDeltasMP[GP(i,j,k,2)] = GlobalMeshW[LW(i,j,k+1,2)] - GlobalMeshW[GW(i,j,k,2)]; // Delta Z
			}
		}
	}

	

	// Staggered U mesh distances
	for(j = 0; j < NY; j++){
		for(k = 0; k < NZ; k++){

			GlobalDeltasMU[GU(0,j,k,0)] = GlobalMeshP[GP(0,j,k,0)] - GlobalMeshU[GU(0,j,k,0)]; // Delta X Parte Izquierda
			GlobalDeltasMU[GU(0,j,k,1)] = GlobalMeshV[GV(0,j+1,k,1)] - GlobalMeshV[GV(0,j,k,1)]; // Delta Y Parte Izquierda
			GlobalDeltasMU[GU(0,j,k,2)] = GlobalMeshW[GW(0,j,k+1,2)] - GlobalMeshW[GW(0,j,k,2)]; // Delta Z Parte Izquierda

			for(i = 1; i < NX; i++){
				GlobalDeltasMU[GU(i,j,k,0)] = GlobalMeshP[GP(i,j,k,0)] - GlobalMeshP[GP(i-1,j,k,0)]; // Delta X
				GlobalDeltasMU[GU(i,j,k,1)] = GlobalMeshV[GV(i,j+1,k,1)] - GlobalMeshV[GV(i,j,k,1)]; // Delta Y
				GlobalDeltasMU[GU(i,j,k,2)] = GlobalMeshW[GW(i,j,k+1,2)] - GlobalMeshW[GW(i,j,k,2)]; // Delta Z
			}

			GlobalDeltasMU[GU(NX,j,k,0)] = GlobalMeshU[GU(NX,j,k,0)] - GlobalMeshP[GP(NX-1,j,k,0)]; //Delta X Parte Derecha
			GlobalDeltasMU[GU(NX,j,k,1)] = GlobalMeshV[GV(NX-1,j+1,k,1)] - GlobalMeshV[GV(NX-1,j,k,1)]; //Delta Y Parte Derecha
			GlobalDeltasMU[GU(NX,j,k,2)] = GlobalMeshW[GW(NX-1,j,k+1,2)] - GlobalMeshW[GW(NX-1,j,k,2)]; //Delta Z Parte Derecha

		}
	}

	// Staggered V mesh distances
	for(i = 0; i < NX; i++){
		for(k = 0; k < NZ; k++){

			GlobalDeltasMV[GV(i,0,k,0)] = GlobalMeshU[GU(i+1,0,k,0)] - GlobalMeshU[GU(i,0,k,0)]; // Delta X Parte Abajo
			GlobalDeltasMV[GV(i,0,k,1)] = GlobalMeshP[GP(i,0,k,1)] - GlobalMeshV[GV(i,0,k,1)]; // Delta Y Parte Abajo
			GlobalDeltasMV[GV(i,0,k,2)] = GlobalMeshW[GW(i,0,k+1,2)] - GlobalMeshW[GW(i,0,k,2)]; // Delta Z Parte Abajo

			for(j = 1; j < NY; j++){
				GlobalDeltasMV[GV(i,j,k,0)] = GlobalMeshU[GU(i+1,j,k,0)] - GlobalMeshU[GU(i,j,k,0)]; //Delta X
				GlobalDeltasMV[GV(i,j,k,1)] = GlobalMeshP[GP(i,j,k,1)] - GlobalMeshP[GP(i,j-1,k,1)]; //Delta Y
				GlobalDeltasMV[GV(i,j,k,2)] = GlobalMeshW[GW(i,j,k+1,2)] - GlobalMeshW[GW(i,j,k,2)]; //Delta Z
			}

			GlobalDeltasMV[GV(i,NY,k,0)] = GlobalMeshU[GU(i+1,NY-1,k,0)] - GlobalMeshU[GU(i,NY-1,k,0)]; // Delta X Parte Arriba
			GlobalDeltasMV[GV(i,NY,k,1)] = GlobalMeshV[GV(i,NY,k,1)] - GlobalMeshP[GP(i,NY-1,k,1)]; // Delta Y Parte Arriba
			GlobalDeltasMV[GV(i,NY,k,2)] = GlobalMeshW[GW(i,NY-1,k+1,2)] - GlobalMeshW[GW(i,NY-1,k,2)]; // Delta Z Parte Arriba

		}
	}

	// Staggered W mesh distances
	for(i = 0; i < NX; i++){
		for(j = 0; j < NY; j++){

			GlobalDeltasMW[GW(i,j,0,0)] = GlobalMeshU[GU(i+1,j,0,0)] - GlobalMeshU[GU(i,j,0,0)]; // Delta X Parte Here
			GlobalDeltasMW[GW(i,j,0,1)] = GlobalMeshV[GV(i,j+1,0,1)] - GlobalMeshV[GV(i,j,0,1)]; // Delta Y Parte Here
			GlobalDeltasMW[GW(i,j,0,2)] = GlobalMeshP[GP(i,j,0,2)] - GlobalMeshW[GW(i,j,0,2)]; // Delta Z Parte Here

			for(k = 1; k < NZ; k++){
				GlobalDeltasMW[GW(i,j,k,0)] = GlobalMeshU[GU(i+1,j,k,0)] - GlobalMeshU[GU(i,j,k,0)]; //Delta X
				GlobalDeltasMW[GW(i,j,k,1)] = GlobalMeshV[GV(i,j+1,k,1)] - GlobalMeshV[GV(i,j,k,1)]; //Delta Y
				GlobalDeltasMW[GW(i,j,k,2)] = GlobalMeshP[GP(i,j,k,2)] - GlobalMeshP[GP(i,j,k-1,2)]; //Delta Z
			}

			GlobalDeltasMW[GW(i,j,NZ,0)] = GlobalMeshU[GU(i+1,j,NZ-1,0)] - GlobalMeshU[GU(i,j,NZ-1,0)]; // Delta X Parte There
			GlobalDeltasMW[GW(i,j,NZ,1)] = GlobalMeshV[GV(i,j+1,NZ-1,1)] - GlobalMeshV[GV(i,j,NZ-1,1)]; // Delta Y Parte There
			GlobalDeltasMW[GW(i,j,NZ,2)] = GlobalMeshW[GW(i,j,NZ,2)] - GlobalMeshP[GP(i,j,NZ-1,2)]; // Delta Z Parte There

		}
	}

}

//Ejecutar todos los procesos del mallador
void Mesher::ExecuteMesher(Memory M1){
char Prob[Problema.length()];
char MeshName[Problema.length() + 5];


	Allocate_MesherMemory(M1); // Memory Allocation
	Get_LocalMeshes(); // Creación de todas las mallas
	Get_LocalColumnsNY(); // Seteo del numero de nodos por columna (Local)
	Get_Deltas(); // Cálculo de las distancias entre nodos en cada una de las matrices
	Get_Surfaces(); // Cálculo de las superficies de cada uno de los volúmenes de control
	Get_Volumes(); // Cálculo de los volúmenes de control de cada volúmen
	
	if(Rango == 0){	
		Get_GlobalMesh(); // Creacion de la malla global
		Get_GlobalColumnsNY(); // Seteo del numero de nodos por columna (Global)
		Get_TotalNodes(); // Calculation of global mesh total nodes
		Get_GlobalDeltas(); // Calculo de las distancias entre nodos en cada matriz (Global)

		strcpy(Prob, Problema.c_str());
		sprintf(MeshName, "%s_Mesh", Prob);
		Get_MeshVTK("Meshes/", MeshName); // Importing the mesh to VTK file

		cout<<"TotalNodes: "<<TotalNodesP<<endl;
		cout<<"Mesh created."<<endl;

	}

}

