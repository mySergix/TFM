//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR MESH IMPORT TO PARAVIEW                               //
//------------------------------------------------------------------------------------------------//

// Function to write a VTK file to import the mesh to Paraview
void Mesher::Get_MeshVTK(string Carpeta, string NombreFile){
int i, j, k;

	ofstream file;
    stringstream InitialName;
    string FinalName;

	InitialName<<"../ParaviewResults/"<<Carpeta<<NombreFile<<".vtk";

	FinalName = InitialName.str();
    file.open(FinalName.c_str());

    file<<"# vtk DataFile Version 2.0"<<endl;
    file<<"GlobalMesh"<<endl;
    file<<"ASCII"<<endl;
    file<<endl;
    file<<"DATASET STRUCTURED_GRID"<<endl;
    file<<"DIMENSIONS"<<"   "<<NX<<"   "<<NY<<"   "<<NZ<<endl;
    file<<endl;
    file<<"POINTS"<<"   "<<NX*NY*NZ<<"   "<<"double"<<endl;
	
	for(k = 0; k < NZ; k++){
		for(j = 0; j < NY; j++){
			for(i = 0; i < NX; i++){	
				file<<GlobalMeshP[GP(i,j,k,0)]<<"   "<<GlobalMeshP[GP(i,j,k,1)]<<"   "<<GlobalMeshP[GP(i,j,k,2)]<<endl;
			}
		}
	}
        
    file<<endl;
	file<<"POINT_DATA"<<"   "<<NX*NY*NZ<<endl;
    file<<"SCALARS "<<"GlobalMesh"<<" double"<<endl;
    file<<"LOOKUP_TABLE"<<"   "<<"GlobalMesh"<<endl;
    file<<endl;
	for(k = 0; k < NZ; k++){
		for(j = 0; j < NY; j++){
			for(i = 0; i < NX; i++){	
				file<<0.0<<" ";
			}
		}
	}

    file.close();

}
