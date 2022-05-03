//------------------------------------------------------------------------------------------------//
//              CPP FILE FOR NON PREMIXED CASE TEMPERATURE BOUNDARY CONDITIONS SETUP              //
//------------------------------------------------------------------------------------------------//

// Function to set all the initial Velocity Boundary Conditions of the simulation
void Solver::Get_StaticBoundaryConditions_Temperatures(Mesher MESH){
int i, j, k;

    // Bottom
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
            if (i < NX_1){
                T.Bottom[BOTTOM(i,0,k)] = T_FlowInlet;
            }
            else if (i >= NX_1 + NX_2){
                T.Bottom[BOTTOM(i,0,k)] = T_FlowInlet;
            }
            else{
                T.Bottom[BOTTOM(i,0,k)] = Twalls;
            }
        }
    }

    // Internal Left
    if (Int_Left){
        i = NX_1 - 1;
		for(j = NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < NY_ColumnMP[i + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			for(k = 0; k < NZ; k++){	
				T.I_Left[ILEFT(0,j,k)] = Twalls;
			}
		}
    }

    // Internal Right
    if (Int_Right){
		i = NX_1 + NX_2 + 1;
		for(j = NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < NY_ColumnMP[i + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			for(k = 0; k < NZ; k++){	
				T.I_Right[IRIGHT(0,j,k)] = Twalls;
			}
		}
	}

}

// Function to update the temperatures boundary conditions
void Solver::Get_UpdateBoundaryConditions_Temperatures(double *T_Matrix){
int i, j, k;

    // Left (Symmetry)
    if (Rango == 0){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                T.Left[LEFT(0,j,k)] = T_Matrix[LP(0,j,k,0)];
            }
        }
    }

    // Right (Symmetry)
    if (Rango == Procesos - 1){ 
        i = NX - 1;
        for(j = NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
            for (k = 0; k < NZ; k++){
                T.Right[RIGHT(NX,j,k)] = T_Matrix[LP(NX-1,j,k,0)];
            }
        }
    } 

    // Top (Null Gradient)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
            T.Top[TOP(i,NY,k)] = T_Matrix[LP(i,NY-1,k,0)];
        }
    }

    // Here (Null Gradient)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            T.Here[HERE(i,j,0)] = T_Matrix[LP(i,j,0,0)];
        }
    }

    // There (Null Gradient)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            T.There[THERE(i,j,NZ)] = T_Matrix[LP(i,j,NZ-1,0)];
        }
    }

}

// Function to set all the corresponding static temperatures halos (no changing)
void Solver::Get_StaticHalos_Temperatures(double *T_Matrix){
int i, j, k;

    // Bottom
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = - Halo; j < 0; j++){
            for (k = 0; k < NZ; k++){
                T_Matrix[LP(i,j,k,0)] = T.Bottom[BOTTOM(i,0,k)];
            }
        }   
    }

    // Internal Left
    if (Int_Left){
        
        for (i = NX_1; i < NX_1 + HP; i++){
            for(j = NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < NY_ColumnMP[i + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			    for(k = 0; k < NZ; k++){	
                    T_Matrix[LP(i,j,k,0)] = T.I_Left[ILEFT(0,j,k)];
			    }
		    }
        }
    }

    // Internal Right
    if (Int_Right){
        for (i = NX_1 + NX_2 - Halo; i < NX_1 + NX_2; i++){
            for(j = NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < NY_ColumnMP[i + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			    for(k = 0; k < NZ; k++){	
                    T_Matrix[LP(i,j,k,0)] = T.I_Right[IRIGHT(0,j,k)];
			    }
		    }
        }
	}
 
}

// Function to update the temperature halos
void Solver::Get_UpdateHalos_Temperatures(double *T_Matrix){
int i, j, k;

    // Left (Symmetry)
    if (Rango == 0){
        for (i = - Halo; i < Ix[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    T_Matrix[LP(i,j,k,0)] = T.Left[LEFT(0,j,k)];   
                }
            }
        }  
    }
   
    // Right (Symmetry)
    if (Rango == Procesos - 1){
        for (i = Fx[Rango]; i < Fx[Rango] + Halo; i++){
            for(j = NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
                for (k = 0; k < NZ; k++){
                    T_Matrix[LP(i,j,k,0)] = T.Right[RIGHT(NX,j,k)];
                }
            }
        }    
    }

    // Top
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = NY; j < NY + Halo; j++){
            for (k = 0; k < NZ; k++){
                T_Matrix[LP(i,j,k,0)] = T.Top[TOP(i,NY,k)];
            }
        }  
    }

    // Here
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for(k = - Halo; k < 0; k++){
                T_Matrix[LP(i,j,k,0)] = T.Here[HERE(i,j,0)];
            }  
        }
    }

    // There
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = NZ; k < NZ + Halo; k++){
                T_Matrix[LP(i,j,k,0)] = T.There[THERE(i,j,NZ)];
            }
        }
    }

}