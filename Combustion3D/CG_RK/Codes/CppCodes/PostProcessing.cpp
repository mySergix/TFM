//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR PARALLEL PROGRAMMING CLASS                            //
//------------------------------------------------------------------------------------------------//

using namespace std;

#include "../HeaderCodes/Memory.h"
#include "../HeaderCodes/ReadData.h"
#include "../HeaderCodes/Parallel.h"
#include "../HeaderCodes/Mesher.h"
#include "../HeaderCodes/PostProcessing.h"

PostProcessing::PostProcessing(Memory M1, ReadData R1, Mesher MESH, Parallel P1){
	
    Problema = MESH.Problema;

    NX = MESH.NX;
    NY = MESH.NY;
    NZ = MESH.NZ;

	Halo = 2;
	HP = 2;

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

// Files of the class
#include "Matrix_Index.cpp"

// Function to write a .VTK to see scalar fields results in Paraview
void PostProcessing::VTK_GlobalScalar3D(string Carpeta, string Variable, string NombreFile, Mesher MESH, double *PropertyMatrix){
int i, j, k;

	ofstream file;
    stringstream InitialName;
    string FinalName;

	InitialName<<"../ParaviewResults/"<<Carpeta<<NombreFile<<".vtk";

	FinalName = InitialName.str();
    file.open(FinalName.c_str());

    file<<"# vtk DataFile Version 2.0"<<endl;
    file<<Variable<<endl;
    file<<"ASCII"<<endl;
    file<<endl;
    file<<"DATASET STRUCTURED_GRID"<<endl;
    file<<"DIMENSIONS"<<"   "<<(NX + 1)<<"   "<<(NY + 1)<<"   "<<(NZ + 1)<<endl;
    file<<endl;
    file<<"POINTS"<<"   "<<(NX + 1)*(NY + 1)*(NZ + 1)<<"   "<<"double"<<endl;
	
	for(k = 0; k < NZ + 1; k++){
		for(j = 0; j < NY + 1; j++){
			for(i = 0; i < NX + 1; i++){
				file<<MESH.GlobalMeshU[GU(i,j,k,0)]<<"   "<<MESH.GlobalMeshV[GV(i,j,k,1)]<<"   "<<MESH.GlobalMeshW[GW(i,j,k,2)]<<endl;
			}
		}
	}
        
    file<<endl;
	file<<"POINT_DATA"<<"   "<<(NX + 1)*(NY + 1)*(NZ + 1)<<endl;
    file<<"SCALARS "<<Variable<<" double"<<endl;
    file<<"LOOKUP_TABLE"<<"   "<<Variable<<endl;
    file<<endl;
	for(k = 0; k < NZ + 1; k++){
		for(j = 0; j < NY + 1; j++){
			for(i = 0; i < NX + 1; i++){
				file<<(1.0 / 6.0) * (PropertyMatrix[GP(i-1,j,k,0)] + PropertyMatrix[GP(i,j-1,k,0)] + PropertyMatrix[GP(i,j,k-1,0)] + PropertyMatrix[GP(i,j,k,0)] + PropertyMatrix[GP(i,j,k,0)] + PropertyMatrix[GP(i,j,k,0)])<<" ";
			}
		}
	}

    file.close();

}

// Function to write a .VTK to see vectorial fields results in Paraview
void PostProcessing::VTK_GlobalVectorial3D(string Carpeta, string Variable, string NombreFile, Mesher MESH, double *Field1, double *Field2, double *Field3){
int i, j, k;

	ofstream file;
    stringstream InitialName;
    string FinalName;

	InitialName<<"../ParaviewResults/"<<Carpeta<<NombreFile<<".vtk";

	FinalName = InitialName.str();
    file.open(FinalName.c_str());

    file<<"# vtk DataFile Version 2.0"<<endl;
    file<<Variable<<endl;
    file<<"ASCII"<<endl;
    file<<endl;
    file<<"DATASET STRUCTURED_GRID"<<endl;
    file<<"DIMENSIONS"<<"   "<<(NX + 1)<<"   "<<(NY + 1)<<"   "<<(NZ + 1)<<endl;
    file<<endl;
    file<<"POINTS"<<"   "<<(NX + 1)*(NY + 1)*(NZ + 1)<<"   "<<"double"<<endl;
	
	for(k = 0; k < NZ + 1; k++){
		for(j = 0; j < NY + 1; j++){
			for(i = 0; i < NX + 1; i++){
				file<<MESH.GlobalMeshU[GU(i,j,k,0)]<<"   "<<MESH.GlobalMeshV[GV(i,j,k,1)]<<"   "<<MESH.GlobalMeshW[GW(i,j,k,2)]<<endl;
			}
		}
	}
        
    file<<endl;
    file<<"POINT_DATA"<<"   "<<(NX + 1)*(NY + 1)*(NZ + 1)<<endl;
    file<<"VECTORS "<<Variable<<" double"<<endl;
    file<<endl;

	for(k = 0; k < NZ + 1; k++){
		for(j = 0; j < NY + 1; j++){
			for(i = 0; i < NX + 1; i++){
				file<<Field1[GU(i,j,k,0)]<<" "<<Field2[GV(i,j,k,0)]<<" "<<Field3[GW(i,j,k,0)]<<endl;
			}
		}
	}

    file.close();

}

// Function to delete all the memory allocated for the PostProcessing
void PostProcessing::Delete_PostProcessingMemory(){
	
	delete [] Ix;
    delete [] Fx;

}