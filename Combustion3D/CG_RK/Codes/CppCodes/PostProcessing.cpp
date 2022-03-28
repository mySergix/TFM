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
	
    Problema = R1.ProblemNumericalData[0];

	//Datos Geométricos del problma
	Xdominio = R1.GeometryData[0];
	Ydominio = R1.GeometryData[1];
	Zdominio = R1.GeometryData[2];

    NX = MESH.NX;
    NY = MESH.NY;
    NZ = MESH.NZ;

	Halo = 2;
	HP = 2;

	Reynolds = R1.ProblemPhysicalData[2];
	Rayleigh = R1.ProblemPhysicalData[3];
	Prandtl = R1.ProblemPhysicalData[5];
	Tleft = R1.ProblemPhysicalData[9];
	Tright = R1.ProblemPhysicalData[10];
	
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

	Mass_Flow = M1.AllocateDouble(NX, 1, 1, 1);
	Bulk_Temperature = M1.AllocateDouble(NX, 1, 1, 1);
	T_Gradient = M1.AllocateDouble(NX, 1, 1, 1);
	Y_Gradient = M1.AllocateDouble(NX, 1, 1, 1);
	Bulk_Fraction = M1.AllocateDouble(NX + 2*HP, 1, NZ + 2*HP, 1);

	Channel_Nusselt = M1.AllocateDouble(NX, 1, 1, 1);
	Channel_Sherwood = M1.AllocateDouble(NX, 1, 1, 1);

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
    file<<"DIMENSIONS"<<"   "<<NX<<"   "<<NY<<"   "<<NZ<<endl;
    file<<endl;
    file<<"POINTS"<<"   "<<NX*NY*NZ<<"   "<<"double"<<endl;
	
	for(k = 0; k < NZ; k++){
		for(j = 0; j < NY; j++){
			for(i = 0; i < NX; i++){
				file<<MESH.GlobalMeshP[GP(i,j,k,0)]<<"   "<<MESH.GlobalMeshP[GP(i,j,k,1)]<<"   "<<MESH.GlobalMeshP[GP(i,j,k,2)]<<endl;
			}
		}
	}
        
    file<<endl;
	file<<"POINT_DATA"<<"   "<<NX*NY*NZ<<endl;
    file<<"SCALARS "<<Variable<<" double"<<endl;
    file<<"LOOKUP_TABLE"<<"   "<<Variable<<endl;
    file<<endl;
	for(k = 0; k < NZ; k++){
		for(j = 0; j < NY; j++){
			for(i = 0; i < NX; i++){
				file<<PropertyMatrix[GP(i,j,k,0)]<<" ";
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
    file<<"DIMENSIONS"<<"   "<<NX<<"   "<<NY<<"   "<<NZ<<endl;
    file<<endl;
    file<<"POINTS"<<"   "<<NX*NY*NZ<<"   "<<"double"<<endl;
	
	for(k = 0; k < NZ; k++){
		for(j = 0; j < NY; j++){
			for(i = 0; i < NX; i++){
				file<<MESH.GlobalMeshP[GP(i,j,k,0)]<<"   "<<MESH.GlobalMeshP[GP(i,j,k,1)]<<"   "<<MESH.GlobalMeshP[GP(i,j,k,2)]<<endl;
			}
		}
	}
        
    file<<endl;
    file<<"POINT_DATA"<<"   "<<NX*NY*NZ<<endl;
    file<<"VECTORS "<<Variable<<" double"<<endl;
    file<<endl;

	for(k = 0; k < NZ; k++){
		for(j = 0; j < NY; j++){	
			for(i = 0; i < NX; i++){
				file<<0.50*(Field1[GU(i,j,k,0)] + Field1[GU(i+1,j,k,0)])<<" "<<0.50*(Field2[GV(i,j,k,0)] + Field2[GV(i,j+1,k,0)])<<" "<<0.50*(Field3[GW(i,j,k,0)] + Field3[GW(i,j,k+1,0)])<<endl;
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

// Function to calculate the Nusselt number along the moist air channel
void PostProcessing::Get_NusseltNumber(Mesher MESH, double *T_Matrix, double *U_Matrix, double *Tbottom, double Twater, double Dh, double Qs){
int i, j, k;
	
	// Calculation of the bulk temperature
	for (i = 0; i < NX; i++){
		Mass_Flow[i] = 0.0;
	}

	// Calculation of mass flow
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){
			for (k = 0; k < NZ; k++){
				Mass_Flow[i] += 0.50 * (U_Matrix[GU(i,j,k,0)] + U_Matrix[GU(i+1,j,k,0)]) * MESH.GlobalSupMP[GP(i,j,k,0)];
			}
		}
	}

	// Calculation of the bulk temperature
	for (i = 0; i < NX; i++){
			Bulk_Temperature[i] = 0.0;
	}

	for (i = 0; i < NX; i++){
		for (k = 0; k < NZ; k++){
			for (j = 0; j < NY; j++){
				Bulk_Temperature[i] += (T_Matrix[GP(i,j,k,0)] * 0.50 * (U_Matrix[GU(i,j,k,0)] + U_Matrix[GU(i+1,j,k,0)])) * MESH.GlobalSupMP[GP(i,j,k,0)];
			}
		}
		Bulk_Temperature[i] = Bulk_Temperature[i] / Mass_Flow[i];
	}

	// Calculation of Temperature Gradient
	for (i = 0; i < NX; i++){
		T_Gradient[i] = 0.0;
	}
	
	for (i = 0; i < NX; i++){
		for (k = 0; k < NZ; k++){
			T_Gradient[i] += ((T_Matrix[GP(i,0,k,0)] - Tbottom[GBOTTOM(i,0,k)]) / MESH.GlobalDeltasMV[GV(i,0,k,1)]) * MESH.GlobalDeltasMP[GP(i,0,k,2)];
		}
		T_Gradient[i] = T_Gradient[i] / MESH.Zdominio;
	}

	// Calculation of the Channel Nusselt Number
	for (i = 0; i < NX; i++){
		Channel_Nusselt[i] = 0.0;
	}

	for (i = 0; i < NX; i++){
		Channel_Nusselt[i] = - Dh * T_Gradient[i] / (Tbottom[GBOTTOM(i,0,NZ/2)] - Bulk_Temperature[i]);	
	}
	
	FILE *fp3;
		fp3 = fopen("../NumericalResults/Channel_Nusselt.txt","w");
			fprintf(fp3,"Coordinate X \t Nusselt \n");
			for (i = 0; i < NX; i++){
				fprintf(fp3, "%f \t %f \n", MESH.GlobalMeshP[GP(i,0,0,0)] / (Reynolds * Prandtl * Dh), Channel_Nusselt[i]);
			}				
		fclose(fp3);

}

// Function to calculate the Sherwood number along the moist air channel
void PostProcessing::Get_SherwoodNumber(Mesher MESH, double *Y_Matrix, double *U_Matrix, double *Ybottom, double Dh){
int i, j, k;

	// Calculation of the bulk temperature
	for (i = 0; i < NX; i++){
		Mass_Flow[i] = 0.0;
	}

	// Calculation of mass flow
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){
			for (k = 0; k < NZ; k++){
				Mass_Flow[i] += 0.50 * (U_Matrix[GU(i,j,k,0)] + U_Matrix[GU(i+1,j,k,0)]) * MESH.GlobalSupMP[GP(i,j,k,0)];
			}
		}
	}
	
	// Calculation of the bulk mass fraction
	for (i = 0; i < NX; i++){
		Bulk_Fraction[i] = 0.0;
	}
	
	for (i = 0; i < NX; i++){
		for (k = 0; k < NZ; k++){
			for (j = 0; j < NY; j++){
				Bulk_Fraction[i] += (Y_Matrix[GP(i,j,k,0)] * 0.50 * (U_Matrix[GU(i,j,k,0)] + U_Matrix[GU(i+1,j,k,0)])) * MESH.GlobalSupMP[GP(i,j,k,0)];
			}
		}
		Bulk_Fraction[i] = Bulk_Fraction[i] / Mass_Flow[i];
	}

	// Calculation of Mass Fraction Gradient
	for (i = 0; i < NX; i++){
		Y_Gradient[i] = 0.0;
	}
	
	for (i = 0; i < NX; i++){
		for (k = 0; k < NZ; k++){
			Y_Gradient[i] += ((Y_Matrix[GP(i,0,k,0)] - Ybottom[GBOTTOM(i,0,k)]) / MESH.GlobalDeltasMV[GV(i,0,k,1)]) * MESH.GlobalDeltasMP[GP(i,0,k,2)];
		}
		Y_Gradient[i] = Y_Gradient[i] / MESH.Zdominio;
	}

	// Calculation of the Channel Nusselt Number
	for (i = 0; i < NX; i++){
		Channel_Sherwood[i] = 0.0;
	}

	for (i = 0; i < NX; i++){
		Channel_Sherwood[i] = - Dh * Y_Gradient[i] / (Ybottom[GBOTTOM(i,0,NZ/2)] - Bulk_Fraction[i]);	
	}

	FILE *fp4;
		fp4 = fopen("../NumericalResults/Channel_Sherwood.txt","w");
			fprintf(fp4,"Coordinate X \t Sherwood \n");
			for (i = 0; i < NX; i++){
				fprintf(fp4, "%f \t %f \n", MESH.GlobalMeshP[GP(i,0,0,0)] / (Reynolds * Prandtl * Dh), Channel_Sherwood[i]);
			}				
		fclose(fp4);

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