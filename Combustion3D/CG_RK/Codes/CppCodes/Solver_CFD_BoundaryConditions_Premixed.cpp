//------------------------------------------------------------------------------------------------//
//                      CPP FILE FOR PREMIXED CASE BOUNDARY CONDITIONS SETUP                      //
//------------------------------------------------------------------------------------------------//

// Function to set all the initial Boundary Conditions of the simulation
void Solver::Get_StaticBoundaryConditions_Velocities(Mesher MESH){
int i, j, k;

    // Bottom
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){

        // Inlet
        if (i < NX_1){
            for (k = 0; k < NZ; k++){
                U.Bottom[BOTTOM(i,0,k)] = 0.0;
                V.Bottom[BOTTOM(i,0,k)] = Uref * (1.0 - pow(MESH.MP[LP(i,0,k,0)] / X_1, 2.0));
                W.Bottom[BOTTOM(i,0,k)] = 0.0;
            }
        }
        // Bottom Walls
        else{
            for (k = 0; k < NZ; k++){
                U.Bottom[BOTTOM(i,0,k)] = 0.0;
                V.Bottom[BOTTOM(i,0,k)] = 0.0;
                W.Bottom[BOTTOM(i,0,k)] = 0.0;
            }
        }
        
    }

    // Int Left
    if (Int_Left == true){
        i = NX_1 - 1;
		for(j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			for(k = 0; k < NZ; k++){	
				U.I_Left[ILEFT(0,j,k)] = 0.0;
                V.I_Left[ILEFT(0,j,k)] = 0.0;
                W.I_Left[ILEFT(0,j,k)] = 0.0;
			}
		}
    }

    // Int Right
    if (Int_Right == true){
		i = NX_1 + NX_2 + 1;
		for(j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			for(k = 0; k < NZ; k++){	
				U.I_Right[IRIGHT(0,j - MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0],k)] = 0.0;
                V.I_Right[IRIGHT(0,j - MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0],k)] = 0.0;
                W.I_Right[IRIGHT(0,j - MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0],k)] = 0.0;
			}
		}
	}

    // Right
    if (Rango == Procesos - 1){
        i = NX - 1;
        for (j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
            for (k = 0; k < NZ; k++){
                U.Right[RIGHT(NX,j,k)] = 0.0;
                V.Right[RIGHT(NX,j,k)] = 0.0;
                W.Right[RIGHT(NX,j,k)] = 0.0;
            }
        }
    }   
    
}

// Function to update the velocity boundary conditions
void Solver::Get_UpdateBoundaryConditions_Velocities(Mesher MESH, double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;

    // Left (Symmmetry)
    if (Rango == 0){

        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Left[LEFT(0,j,k)] = U_Matrix[LU(0,j,k,0)];
                V.Left[LEFT(0,j,k)] = V_Matrix[LV(0,j,k,0)];
                W.Left[LEFT(0,j,k)] = W_Matrix[LW(0,j,k,0)];
            }
        }

    }

    // Top
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            U.Top[TOP(i,NY,k)] = U_Matrix[LU(i,NY - 1,k,0)];
            V.Top[TOP(i,NY,k)] = V_Matrix[LV(i,NY - 1,k,0)];
            W.Top[TOP(i,NY,k)] = W_Matrix[LW(i,NY - 1,k,0)];
        }
    }

    // Here
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
            U.Here[HERE(i,j,0)] = U_Matrix[LU(i,j,0,0)];
            V.Here[HERE(i,j,0)] = V_Matrix[LV(i,j,0,0)];
        }
    }

    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
            W.Here[HERE(i,j,0)] = W_Matrix[LW(i,j,2,0)];
        }
    }

    // There
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
            U.There[THERE(i,j,NZ)] = U_Matrix[LU(i,j,NZ-1,0)];
            V.There[THERE(i,j,NZ)] = V_Matrix[LV(i,j,NZ-1,0)];
        }
    }

    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
            W.There[THERE(i,j,NZ)] = W_Matrix[LW(i,j,NZ-2,0)];
        }
    }

}

// Function to update the predictor velocities boundary conditions
void Solver::Get_UpdateBoundaryConditions_PredictorVelocities(Mesher MESH){
int i, j, k;

    // Left (Symmetry)
    if (Rango == 0){
        i = 0;
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                //U.Predictor[LU(0,j,k,0)] = U.Pres[LU(0,j,k,0)];
            }
        }
    }

    // Right
    if (Rango == Procesos - 1){
        i = NX;
        for (j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1]; j++){
            for (k = 0; k < NZ; k++){
                //U.Predictor[LU(NX,j,k,0)] = U.Pres[LU(NX,j,k,0)];
            }
        }
    }

    // Bottom
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
			V.Predictor[LV(i,MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0],k,0)] = V.Bottom[BOTTOM(i,0,k)];
        }
    }

    // Top
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
			V.Predictor[LV(i,NY,k,0)] = V.Pres[LV(i,NY,k,0)];
        }
    }

    // Here (Null Gradient)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMW[i + Halo - Ix[Rango]][1]; j++){
            //W.Predictor[LW(i,j,0,0)] = W.Predictor[LW(i,j,1,0)];
        }
    }

    // There (Null Gradient)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMW[i + Halo - Ix[Rango]][1]; j++){
            //W.Predictor[LW(i,j,NZ,0)] = W.Predictor[LW(i,j,NZ-1,0)];
        }
    }

}

// Function to set all the corresponding static velocity halos (no changing)
void Solver::Get_StaticHalos_Velocity(Mesher MESH, double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;

    // Bottom
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0] - Halo; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j++){
            for (k = 0; k < NZ; k++){
                    U_Matrix[LU(i,j,k,0)] = U.Bottom[BOTTOM(i,0,k)];
                    W_Matrix[LW(i,j,k,0)] = W.Bottom[BOTTOM(i,0,k)];
            }
        }
    }
    
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] - Halo; j <= MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0]; j++){
            for (k = 0; k < NZ; k++){
                    V_Matrix[LV(i,j,k,0)] = V.Bottom[BOTTOM(i,0,k)];
            }
        }
    }

    // Int Left
    if (Int_Left){
        
        for (i = NX_1; i < NX_1 + Halo; i++){
            for(j = 0; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			    for(k = 0; k < NZ; k++){	
                    V_Matrix[LV(i,j,k,0)] = V.I_Left[ILEFT(0,j,k)];
                    W_Matrix[LW(i,j,k,0)] = W.I_Left[ILEFT(0,j,k)];
			    }
		    }
        }
		
        for (i = NX_1; i < NX_1 + Halo + 1; i++){
            for(j = 0; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			    for(k = 0; k < NZ; k++){	
				    U_Matrix[LU(i,j,k,0)] = U.I_Left[ILEFT(0,j,k)];
			    }
		    }
        }

    }

    // Int Right
    if (Int_Right){
        for (i = NX_1 + NX_2 - Halo; i < NX_1 + NX_2; i++){
            for(j = MESH.NY_ColumnMP[NX_1 + NX_2 + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[NX_1 + NX_2 + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			    for(k = 0; k < NZ; k++){	
                    V_Matrix[LV(i,j,k,0)] = V.I_Right[IRIGHT(0,j - MESH.NY_ColumnMP[NX_1 + NX_2 + Halo - Ix[Rango]][0],k)];
                    W_Matrix[LW(i,j,k,0)] = W.I_Right[IRIGHT(0,j - MESH.NY_ColumnMP[NX_1 + NX_2 + Halo - Ix[Rango]][0],k)];
			    }
		    }
        }
		
		for (i = NX_1 + NX_2 - Halo; i <= NX_1 + NX_2; i++){
            for(j = MESH.NY_ColumnMP[NX_1 + NX_2 + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[NX_1 + NX_2 + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			    for(k = 0; k < NZ; k++){	
				    U_Matrix[LU(i,j,k,0)] = U.I_Right[IRIGHT(0,j - MESH.NY_ColumnMP[NX_1 + NX_2 + Halo - Ix[Rango]][0],k)];
			    }
		    }
        }
	}

    // Right
    if (Rango == Procesos - 1){
        for (i = Fx[Rango]; i < Fx[Rango] + Halo; i++){
            for (j = MESH.NY_ColumnMP[NX - 1 + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[NX - 1 + Halo - Ix[Rango]][1]; j++){
                for (k = 0; k < NZ; k++){
                    V_Matrix[LV(i,j,k,0)] = V.Right[RIGHT(NX,j,k)];
                    W_Matrix[LW(i,j,k,0)] = W.Right[RIGHT(NX,j,k)];
                }
            }
        }

        for (i = Fx[Rango]; i < Fx[Rango] + Halo + 1; i++){
            for (j = MESH.NY_ColumnMU[NX + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMU[NX + Halo - Ix[Rango]][1]; j++){
                for (k = 0; k < NZ; k++){
                    U_Matrix[LU(i,j,k,0)] = U.Right[RIGHT(NX,j,k)];
                }
            }
        }

    }     

}

// Function to update the boundary conditions
void Solver::Get_UpdateHalos_Velocity(Mesher MESH, double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;

    // Left
    if (Rango == 0){
        for (i = - Halo; i < 0; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){     
                    V_Matrix[LV(i,j,k,0)] = V.Left[LEFT(0,j,k)];
                    W_Matrix[LW(i,j,k,0)] = W.Left[LEFT(0,j,k)];
                }
            }
        }

        for (i = - Halo; i <= 0; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){     
                    U_Matrix[LU(i,j,k,0)] = U_Matrix[LU(-i,j,k,0)];
                }
            }
        }

    }

    // Top
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = NY; j < NY + Halo; j++){
            for(k = 0; k < NZ; k++){
                U_Matrix[LU(i,j,k,0)] = U.Top[TOP(i,NY,k)];
                W_Matrix[LW(i,j,k,0)] = W.Top[TOP(i,NY,k)];
            }
        }
    }

    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = NY; j < NY + Halo + 1; j++){
            for(k = 0; k < NZ; k++){
                V_Matrix[LV(i,j,k,0)] = V.Top[TOP(i,NY,k)];
            }
        }
    }

    // Here
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
            for(k = - Halo; k < 0; k++){
                U_Matrix[LU(i,j,k,0)] = U.Here[HERE(i,j,0)];
                V_Matrix[LV(i,j,k,0)] = V.Here[HERE(i,j,0)];
            }  
        }
    }

    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
            for(k = - Halo; k <= 1; k++){
                W_Matrix[LW(i,j,k,0)] = W.Here[HERE(i,j,0)];
            }  
        }
    }

    // There
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
            for (k = NZ; k < NZ + Halo; k++){
                U_Matrix[LU(i,j,k,0)] = U.There[THERE(i,j,NZ)];
                V_Matrix[LV(i,j,k,0)] = V.There[THERE(i,j,NZ)];
            } 
        }
    }

    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
            for (k = NZ-1; k < NZ + Halo + 1; k++){
                W_Matrix[LW(i,j,k,0)] = W.There[THERE(i,j,NZ)];
            } 
        }
    }

}