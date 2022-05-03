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

	Rayleigh = ;
	Tleft = ;
	Tright = ;
	
	//Datos necesarios para computación paralela
	Rango = P1.Rango;
	Procesos = P1.Procesos;
	Ix = M1.AllocateInt(Procesos, 1, 1, 1);
    Fx = M1.AllocateInt(Procesos, 1, 1, 1);

    for (int i = 0; i < Procesos; i++){
        Ix[i] = P1.Ix[i];
        Fx[i] = P1.Fx[i];
    }

	LocalNusselt = M1.AllocateDouble(NY, 1, 1, 1);
	U_Averaged = M1.AllocateDouble(NY, 1, 1, 1);
	V_Averaged = M1.AllocateDouble(NX, 1, 1, 1);

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
				file<<GlobalMeshU[GU(i,j,k,0)]<<"   "<<GlobalMeshV[GV(i,j,k,1)]<<"   "<<GlobalMeshW[GW(i,j,k,2)]<<endl;
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
				file<<GlobalMeshU[GU(i,j,k,0)]<<"   "<<GlobalMeshV[GV(i,j,k,1)]<<"   "<<GlobalMeshW[GW(i,j,k,2)]<<endl;
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

// Function to calculate the Nusselt number at the walls
void PostProcessing::Get_NusseltResults(Mesher MESH, double *T_matrix){
int i, j, k;
double NusseltMax, NusseltMin, NusseltMedio;
double Umax = 0.0;
int ImaxU = 0;
double Vmax = 0.0;
int ImaxV = 0;
double Sumatorio = 0.0;
int ImaxNu = 0;
int IminNu = 0;

	for (j = 0; j < NY; j++){
		LocalNusselt[j] = 0.0;
	}

	for (j = 0; j <  NY; j++){
		for (k = 0; k < NZ; k++){
			LocalNusselt[j] += ((Tleft - T_matrix[GP(0,j,k,0)]) / MESH.DeltasMU[LU(0,j,k,0)]) * (MESH.DeltasMP[LP(0,j,k,2)] / MESH.Zdominio);
		}
		LocalNusselt[j] = LocalNusselt[j] * (MESH.Xdominio / (Tleft - Tright));
		Sumatorio += LocalNusselt[j] * MESH.DeltasMP[LP(0,j,NZ/2,1)];
	} 

	NusseltMedio = Sumatorio / MESH.Ydominio;

	NusseltMax = 0.0;
	NusseltMin = 1e3;

	for (j = 0; j < NY; j++){
		if (LocalNusselt[j] >= NusseltMax){
			NusseltMax = LocalNusselt[j];
			ImaxNu = j;
		}
	}

	for (j = 0; j < NY; j++){
		if (LocalNusselt[j] <= NusseltMin){
			NusseltMin = LocalNusselt[j];
			IminNu = j;
		}
	}

	FILE *fp1;
		fp1 = fopen("../NumericalResults/NusseltResults.txt","w");
			fprintf(fp1,"Número de Rayleigh: %f \n", Rayleigh);
			fprintf(fp1, "\n");
			fprintf(fp1,"Número Nusselt medio: %f \n", NusseltMedio);
			fprintf(fp1, "\n");
			fprintf(fp1,"Número de Nusselt máximo: %f \n", NusseltMax);
			fprintf(fp1,"Posición Y del número de Nusselt máximo (m): %f \n", MESH.MP[LP(0,ImaxNu,NZ/2,1)]);
			fprintf(fp1, "\n");
			fprintf(fp1,"Número de Nusselt mínimo: %f \n", NusseltMin);
			fprintf(fp1,"Posición Y del número de Nusselt mínimo (m): %f \n", MESH.MP[LP(0,IminNu,NZ/2,1)]);
					
		fclose(fp1);

}

// Function to calculate the average velocity at the mid sections
void PostProcessing::Get_VelocityResults(Mesher MESH, double *U_matrix, double *V_matrix, double K){
int i, j, k;

	// U Velocity Mid Vertical Line (Z - Averaged)
	for (j = 0; j < NY; j++){
		U_Averaged[j] = 0.0;
	}

	for (j = 0; j < NY; j++){
		for (k = 0; k < NZ; k++){
			U_Averaged[j] += U_matrix[GU(NX/2,j,k,0)] * (MESH.DeltasMU[LU(0,j,k,2)] / MESH.Zdominio);
		}
	}

	// V Velocity Mid Horizonal Line (Z - Averaged)
	for (i = 0; i < NX; i++){
		V_Averaged[i] = 0.0;
	}

	for (i = 0; i < NX; i++){
		for (k = 0; k < NZ; k++){
			V_Averaged[i] += V_matrix[GV(i,NY/2,k,0)] * (MESH.GlobalDeltasMV[GV(i,NY/2,k,2)] / MESH.Zdominio);
		}
	}

	if (Problema == 1){
		FILE *fp1;
		fp1 = fopen("../NumericalResults/Driven_U_Results.txt","w");
			fprintf(fp1,"Coordinate \t U Velocity \n");
			for (j = 0; j < NY; j++){
				fprintf(fp1, "%f \t %f \n", MESH.GlobalMeshU[GU(NX/2,j,0,1)], U_Averaged[j]);
			}				
		fclose(fp1);

		FILE *fp2;
		fp2 = fopen("../NumericalResults/Driven_V_Results.txt","w");
			fprintf(fp2,"Coordinate \t V Velocity \n");
			for (i = 0; i < NX; i++){
				fprintf(fp2, "%f \t %f \n", MESH.GlobalMeshV[GV(i,NY/2,0,0)], V_Averaged[i]);
			}				
		fclose(fp2);
	}
	else if (Problema == 2){
		FILE *fp1;
		fp1 = fopen("../NumericalResults/Differentially_U_Results.txt","w");
			fprintf(fp1,"Coordinate \t U Velocity \n");
			for (j = 0; j < NY; j++){
				fprintf(fp1, "%f \t %f \n", MESH.GlobalMeshU[GU(NX/2,j,0,1)], U_Averaged[j] * (MESH.Xdominio / K));
			}				
		fclose(fp1);

		FILE *fp2;
		fp2 = fopen("../NumericalResults/Differentially_V_Results.txt","w");
			fprintf(fp2,"Coordinate \t V Velocity \n");
			for (i = 0; i < NX; i++){
				fprintf(fp2, "%f \t %f \n", MESH.GlobalMeshV[GV(i,NY/2,0,0)], V_Averaged[i] * (MESH.Xdominio / K));
			}				
		fclose(fp2);
	}

}

// Function to delete all the memory allocated for the PostProcessing
void PostProcessing::Delete_PostProcessingMemory(){

	if (Rango == 0){
		delete [] LocalNusselt;
		delete [] U_Averaged;
		delete [] V_Averaged;
	}
	
	delete [] Ix;
    delete [] Fx;

}