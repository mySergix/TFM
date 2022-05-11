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

    NX_1 = MESH.NX_1;
	NX_2 = MESH.NX_2;
	NX_3 = MESH.NX_3;

	NY_1 = MESH.NY_1;
	NY_2 = MESH.NY_2;
	NY_3 = MESH.NY_3;
	NY_4 = MESH.NY_4;

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

// Function to set the global halos of a scalar matrix
void PostProcessing::Get_GlobalScalarHalos(Mesher MESH, double *ScalarMatrix){
int i, j, k;

    // Right
    for (j = 0; j < NY; j++){
        for (k = 0; k < NZ; k++){
            ScalarMatrix[GP(NX,j,k,0)] = ScalarMatrix[GP(NX-1,j,k,0)];
        }
    }

    // Top
    for (i = 0; i < NX; i++){
        for (k = 0; k < NZ; k++){
            ScalarMatrix[GP(i,NY,k,0)] = ScalarMatrix[GP(i,NY-1,k,0)];
        }
    }

    // There
    for (i = 0; i < NX; i++){
        for (j = 0; j < NY; j++){
            ScalarMatrix[GP(i,j,NZ,0)] = ScalarMatrix[GP(i,j,NZ-1,0)];
        }
    }

    // Top - There Corner
    for (i = 0; i < NX; i++){
        ScalarMatrix[GP(i,NY,NZ,0)] = 0.50 * (ScalarMatrix[GP(i,NY-1,NZ,0)] + ScalarMatrix[GP(i,NY,NZ-1,0)]);
    }

    // Right - Top Corner
    for (k = 0; k < NZ; k++){
        ScalarMatrix[GP(NX,NY,k,0)] = 0.50 * (ScalarMatrix[GP(NX-1,NY,k,0)] + ScalarMatrix[GP(NX,NY-1,k,0)]);
    }

    // Right - There Corner
    for (j = 0; j < NY; j++){
        ScalarMatrix[GP(NX,j,NZ,0)] = 0.50 * (ScalarMatrix[GP(NX-1,j,NZ,0)] + ScalarMatrix[GP(NX,j,NZ-1,0)]);
    }

    // Top - There - Right Corner
    ScalarMatrix[GP(NX,NY,NZ,0)] = (1.0 / 3.0) * (ScalarMatrix[GP(NX-1,NY,NZ,0)] + ScalarMatrix[GP(NX,NY,NZ-1,0)] + ScalarMatrix[GP(NX,NY-1,NZ,0)]);

    // Solid Walls

    // Int Left
    for (j = MESH.GlobalNY_ColumnMP[NX_1 - 1][0]; j < MESH.GlobalNY_ColumnMP[NX_1 - 1][1]; j++){
        for (k = 0; k < NZ + 1; k++){
            ScalarMatrix[GP(NX_1,j,k,0)] = ScalarMatrix[GP(NX_1-1,j,k,0)];
        }
    }

    // Int Right
    for (j = MESH.GlobalNY_ColumnMP[NX_1 + NX_2][0]; j < MESH.GlobalNY_ColumnMP[NX_1 + NX_2][1]; j++){
        for (k = 0; k < NZ + 1; k++){
            ScalarMatrix[GP(NX_1 + NX_2 - 1,j,k,0)] = ScalarMatrix[GP(NX_1 + NX_2,j,k,0)];
        }
    }
    
    // Slit Wall
    for (i = NX_1; i < NX_1 + NX_2; i++){
        for (k = 0; k < NZ; k++){
            ScalarMatrix[GP(i,MESH.GlobalNY_ColumnMP[i][0]-1,k,0)] = ScalarMatrix[GP(i,MESH.GlobalNY_ColumnMP[i][0],k,0)];
        }
    }

}

// Function to set the global halos of the vectorial matrix
void PostProcessing::Get_GlobalVectorialHalos(double *U_matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;

    // Velocity U


    // Velocity V

    // Right
    for (j = 0; j < NY + 1; j++){
        for (k = 0; k < NZ + 1; k++){
            V_Matrix[GV(NX,j,k,0)] = V_Matrix[GV(NX-1,j,k,0)];
        }
    }

    // There
    for (i = 0; i < NX + 1; i++){
        for (j = 0; j < NY + 1; j++){
            V_Matrix[GV(i,j,NZ,0)] = V_Matrix[GV(i,j,NZ-1,0)];
        }
    }


}

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
                file<<PropertyMatrix[GP(i,j,k,0)]<<" ";
				//file<<(1.0 / 6.0) * (PropertyMatrix[GP(i-1,j,k,0)] + PropertyMatrix[GP(i,j-1,k,0)] + PropertyMatrix[GP(i,j,k-1,0)] + PropertyMatrix[GP(i,j,k,0)] + PropertyMatrix[GP(i,j,k,0)] + PropertyMatrix[GP(i,j,k,0)])<<" ";
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