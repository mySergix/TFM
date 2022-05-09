//------------------------------------------------------------------------------------------------//
//                   CPP FILE FOR CFD SOLVER CLASS POISSON COEFFICIENTS                           //
//------------------------------------------------------------------------------------------------//

// Function to calculate the Poisson Coefficients
void Solver::Get_PoissonCoefficients(Mesher MESH){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                A.aw[LA(i,j,k,0)] = 0.0;
                A.ae[LA(i,j,k,0)] = 0.0; 

                A.as[LA(i,j,k,0)] = 0.0; 
                A.an[LA(i,j,k,0)] = 0.0; 

                A.ah[LA(i,j,k,0)] = 0.0; 
                A.at[LA(i,j,k,0)] =  0.0;
            }
        }
    }

    // West and East Coefficients (aw, ae)

    if (Problema == "Premixed"){
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
		        for(k = 0; k < NZ; k++){	

                    // Core
                    A.aw[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,0)]/MESH.DeltasMU[LU(i,j,k,0)];
                    A.ae[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,0)]/MESH.DeltasMU[LU(i+1,j,k,0)];

                    // Internal Left
                    if (i == NX_1 - 1 && j < NY_1 + NY_2){
                        A.ae[LA(i,j,k,0)] = 0.0;
                    }
                
                    // Internal Right
                    if (i == NX_1 + NX_2 && j >= NY_1 && j < NY_1 + NY_2){
                        A.aw[LA(i,j,k,0)] = 0.0;
                    }

                }
            }
        }
    }
    else if (Problema == "NonPremixed"){
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = 0; j < NY; j++){
		        for(k = 0; k < NZ; k++){	

                    // Core
                    A.aw[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,0)]/MESH.DeltasMU[LU(i,j,k,0)];
                    A.ae[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,0)]/MESH.DeltasMU[LU(i+1,j,k,0)];

                    // Internal Left
                    if (i == NX_1 - 1 && j < NY_2){
                        A.ae[LA(i,j,k,0)] = 0.0;
                    }
                
                    // Internal Right
                    if (i == NX_1 + NX_2 && j < NY_2){
                        A.aw[LA(i,j,k,0)] = 0.0;
                    }

                }
            }
        }
    }
    
    if (Rango == 0){
        for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
                A.aw[LA(0,j,k,0)] = 0.0;
            }
        }
    }
    else if (Rango == Procesos - 1){
        for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
                A.ae[LA(NX-1,j,k,0)] = 0.0;
            }
        }
    }

    // South and North Coefficients (as, an)
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(k = 0; k < NZ; k++){
            
            // Top
            A.as[LA(i,NY-1,k,0)] = MESH.SupMP[LP(i,NY-1,k,1)]/MESH.DeltasMV[LV(i,NY-1,k,1)];
            A.an[LA(i,NY-1,k,0)] = 0.0;

            // Bottom
            A.as[LA(i,MESH.NY_ColumnMP[i + HP - Ix[Rango]][0],k,0)] = 0.0;
            A.an[LA(i,MESH.NY_ColumnMP[i + HP - Ix[Rango]][0],k,0)] = MESH.SupMP[LP(i,MESH.NY_ColumnMP[i + HP - Ix[Rango]][0],k,1)]/MESH.DeltasMV[LV(i,MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1,k,1)];

			for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){

                // Core
                A.as[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,1)]/MESH.DeltasMV[LV(i,j,k,1)];
                A.an[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,1)]/MESH.DeltasMV[LV(i,j+1,k,1)];

            }
        }
    }

    // Here and There Coefficients (ah, at)
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){

            // Here
            A.ah[LA(i,j,0,0)] = 0.0;//MESH.SupMP[LP(i,j,0,2)] / (MESH.DeltasMW[LW(i,j,NZ,2)] + MESH.DeltasMW[LW(i,j,0,2)]);
			A.at[LA(i,j,0,0)] = MESH.SupMP[LP(i,j,0,2)]/MESH.DeltasMW[LW(i,j,1,2)];

            // There
            A.ah[LA(i,j,NZ-1,0)] = MESH.SupMP[LP(i,j,NZ-1,2)]/MESH.DeltasMW[LW(i,j,NZ-1,2)];
			A.at[LA(i,j,NZ-1,0)] = 0.0;//MESH.SupMP[LP(i,j,NZ,2)] / (MESH.DeltasMW[LW(i,j,0,2)] + MESH.DeltasMW[LW(i,j,NZ,2)]);

            for(k = 1; k < NZ - 1; k++){
                A.ah[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,2)]/MESH.DeltasMW[LW(i,j,k,2)];
				A.at[LA(i,j,k,0)] = MESH.SupMP[LP(i,j,k,2)]/MESH.DeltasMW[LW(i,j,k+1,2)];
            }
        }
    }

	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(k = 0; k < NZ; k++){
			for(j = 0; j < NY; j++){
                A.ap[LA(i,j,k,0)] = A.aw[LA(i,j,k,0)] + A.ae[LA(i,j,k,0)] + A.as[LA(i,j,k,0)] + A.an[LA(i,j,k,0)] + A.ah[LA(i,j,k,0)] + A.at[LA(i,j,k,0)] + 1e-12;
            }
        }
    }

}

// Function to calculate the total number of nonzero elements in Laplacian Matrix
void Solver::Get_NonZero_NumberElements(){
int i, j, k;

    // Number of nonzero elements
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){

                // Aw Coefficient
                if(A.aw[LA(i,j,k,0)] > 0.0){
                    NNZ++;
                }

                // Ae Coefficient
                if(A.ae[LA(i,j,k,0)] > 0.0){
                    NNZ++;
                }

                // As Coefficient
                if(A.as[LA(i,j,k,0)] > 0.0){
                    NNZ++;
                }

                // An Coefficient
                if(A.an[LA(i,j,k,0)] > 0.0){
                    NNZ++;
                }

                // Ah Coefficient
                if(A.ah[LA(i,j,k,0)] > 0.0){
                    NNZ++;
                }

                // At Coefficient
                if(A.at[LA(i,j,k,0)] > 0.0){
                    NNZ++;
                }

                // Ap Coefficient
                if(A.ap[LA(i,j,k,0)] > 0.0){
                    NNZ++;
                }

            }
        }
    }

}

// Function to calculate the global index of the Laplacian Coefficients
int Solver::Get_GlobalIndCoefficient(Mesher MESH, int i, int j, int k){
int Coefficient;

    if (i < NX_1){
        Coefficient = (NY) * (NZ) * (i) + (NZ) * (j) + (k);
	}
	else if (i >= NX_1 && i < NX_1 + NX_2){
        Coefficient = (NY)*(NZ)*(NX_1) + (NY - MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]) * (NZ) * ((i) - NX_1) + (NZ) * ((j) - MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]) + (k);
	}
	else if (i >= NX_1 + NX_2){
		Coefficient = (NY)*(NZ)*(NX_1) + (NY_3 + NY_4)*(NZ)*(NX_2) + (NY - MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]) * (NZ) * ((i) - (NX_1 + NX_2)) + (NZ) * ((j) - MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]) + (k);
	}

    return Coefficient;

}

// Function to calculate the global index of the Laplacian Coefficients
int Solver::Get_LocalIndCoefficient(Mesher MESH, int I, int J, int K){
int LocalCoefficient = 0;
int i;

    for (i = Ix[Rango]; i < I; i++){
        LocalCoefficient += (MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]) * (NZ);
    }

    LocalCoefficient += (J - MESH.NY_ColumnMP[I + HP - Ix[Rango]][0]) * (NZ) + (K);

    return LocalCoefficient;

}

// Function to create the CSR Matrix of the Laplacian (Poisson)
void Solver::Get_CSR_LaplacianMatrix(Mesher MESH){
int i, j, k;
int Val_Ind = 0;
int Row_Count = 0;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
            for (k = 0; k < NZ; k++){
                
                // Check Aw
                if (A.aw[LA(i,j,k,0)] > 0.0){
                    Val_Laplacian[Val_Ind] = A.aw[LA(i,j,k,0)];
                    Col_Ind[Val_Ind] = Get_GlobalIndCoefficient(MESH, i-1, j, k); 
                    Row_Ptr[Row_Count] = Val_Ind;
                    Row_Count++; 
                    Val_Ind++;    
                }

                // Check As
                if (A.as[LA(i,j,k,0)] > 0.0){
                    Val_Laplacian[Val_Ind] = A.as[LA(i,j,k,0)];
                    Col_Ind[Val_Ind] = Get_GlobalIndCoefficient(MESH, i, j-1, k);
                    if (A.aw[LA(i,j,k,0)] == 0.0){
                        Row_Ptr[Row_Count] = Val_Ind;
                        Row_Count++;
                    }
                    Val_Ind++;     
                }

                // Check Ah
                if (A.ah[LA(i,j,k,0)] > 0.0){
                    Val_Laplacian[Val_Ind] = A.ah[LA(i,j,k,0)];
                    Col_Ind[Val_Ind] = Get_GlobalIndCoefficient(MESH, i, j, k-1);
                    if (A.aw[LA(i,j,k,0)] == 0.0 && A.as[LA(i,j,k,0)] == 0.0){
                        Row_Ptr[Row_Count] = Val_Ind;
                        Row_Count++;
                    }
                    Val_Ind++;     
                }

                // Check Ap (would be the first if all previous ones are zero)
                if (A.ap[LA(i,j,k,0)] > 0.0){
                    Val_Laplacian[Val_Ind] = - A.ap[LA(i,j,k,0)];
                    Col_Ind[Val_Ind] = Get_GlobalIndCoefficient(MESH, i, j, k);
                    if (A.aw[LA(i,j,k,0)] == 0.0 && A.as[LA(i,j,k,0)] == 0.0 && A.ah[LA(i,j,k,0)] == 0.0){
                        Row_Ptr[Row_Count] = Val_Ind;
                        Row_Count++;
                    }
                    Val_Ind++;     
                }

                // Check At
                if (A.at[LA(i,j,k,0)] > 0.0){
                    Val_Laplacian[Val_Ind] = A.at[LA(i,j,k,0)];
                    Col_Ind[Val_Ind] = Get_GlobalIndCoefficient(MESH, i, j, k+1);
                    Val_Ind++;     
                }

                // Check An
                if (A.an[LA(i,j,k,0)] > 0.0){
                    Val_Laplacian[Val_Ind] = A.an[LA(i,j,k,0)];
                    Col_Ind[Val_Ind] = Get_GlobalIndCoefficient(MESH, i, j+1, k);
                    Val_Ind++;     
                }

                // Check Ae
                if (A.ae[LA(i,j,k,0)] > 0.0){
                    Val_Laplacian[Val_Ind] = A.ae[LA(i,j,k,0)];
                    Col_Ind[Val_Ind] = Get_GlobalIndCoefficient(MESH, i+1, j, k);
                    Val_Ind++;     
                }

            }
        }
    }

    // Last element of Row_Ptr is the total number of NNZ
    Row_Ptr[Row_Count] = Val_Ind;
    
}

// Function to calculate the index of the RHS vector (bp)
void Solver::Get_RHS_VectorIndex(Mesher MESH){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
            for (k = 0; k < NZ; k++){
                RHS_Ind[Get_LocalIndCoefficient(MESH, i, j, k)] = Get_GlobalIndCoefficient(MESH, i, j, k);
            }
        }
    }

}

// Function to retrieve the solution from PETSc with local Index
void Solver::Get_LocalSolution(Mesher MESH){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
            for (k = 0; k < NZ; k++){
                P.Pres[LP(i,j,k,0)] = X_Sol_Array[Get_LocalIndCoefficient(MESH, i, j, k)];
            }
        }
    }

}