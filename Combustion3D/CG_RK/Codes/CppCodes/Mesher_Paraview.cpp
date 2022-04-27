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
    file<<"SCALARS "<<"GlobalMesh"<<" double"<<endl;
    file<<"LOOKUP_TABLE"<<"   "<<"GlobalMesh"<<endl;
    file<<endl;
	for(k = 0; k < NZ + 1; k++){
		for(j = 0; j < NY + 1; j++){
			for(i = 0; i < NX + 1; i++){	
				file<<0.0<<" ";
			}
		}
	}

    file.close();

}
