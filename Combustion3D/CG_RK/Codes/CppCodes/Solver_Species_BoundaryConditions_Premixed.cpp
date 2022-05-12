//------------------------------------------------------------------------------------------------//
//                 CPP FILE FOR PREMIXED CASE SPECIES BOUNDARY CONDITIONS SETUP                   //
//------------------------------------------------------------------------------------------------//

// Function to set all the initial Species Boundary Conditions of the simulation
void Solver::Get_Species_StaticBoundaryConditions(Mesher MESH){
int i, j, k, n;

    // Inlet Flow (Stoichiometric Conditions)
    if (Ix[Rango] <= NX_1){
        for (i = Ix[Rango]; i < min(Fx[Rango], NX_1); i++){
            for (k = 0; k < NZ; k++){

                // Species 0 - CH4
                Species[0].Bottom[BOTTOM(i,0,k)] = Species[0].Wmolar / (Species[0].Wmolar + 2.0 * Species[1].Wmolar + 2.0 * 3.76 * Species[2].Wmolar);

                // Species 1 - O2
                Species[1].Bottom[BOTTOM(i,0,k)] = Species[1].Wmolar / (Species[0].Wmolar + 2.0 * Species[1].Wmolar + 2.0 * 3.76 * Species[2].Wmolar);

                // Species 2 - N2
                Species[2].Bottom[BOTTOM(i,0,k)] = Species[2].Wmolar / (Species[0].Wmolar + 2.0 * Species[1].Wmolar + 2.0 * 3.76 * Species[2].Wmolar);

                // Species 3 - CO2
                Species[3].Bottom[BOTTOM(i,0,k)] = 0.0;

                // Species 4 - H2O
                Species[4].Bottom[BOTTOM(i,0,k)] = 0.0;

            }
        }
    }
  
}

// Function to update the temperatures boundary conditions
void Solver::Get_Species_UpdateBoundaryConditions(Mesher MESH){
int i, j, k, n;

    // Top (Outlet)
    for (n = 0; n < N_Species; n++){

        for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
            for (k = 0; k < NZ; k++){
                Species[n].Top[TOP(i,NY,k)] = Species[n].Y_Pres[LP(i,NY-1,k,0)];
            }
        }

    }
    
    // Left (Symmetry)
    if (Rango == 0){
        for (n = 0; n < N_Species; n++){

            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    Species[n].Left[LEFT(0,j,k)] = Species[n].Y_Pres[LP(0,j,k,0)];
                }
            }

        }  
    }

    // Bottom (Inlet + Slit + Wall)
    for (n = 0; n < N_Species; n++){

        if (Fx[Rango] >= NX_1){
            for (i = max(Ix[Rango],NX_1); i < Fx[Rango] + 1; i++){
                for (k = 0; k < NZ; k++){
                    Species[n].Bottom[BOTTOM(i,0,k)] = Species[n].Y_Pres[LP(i,MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1,k,0)];
                }
            }
        }
        
    }

    // Right (Burner Wall)
    if (Rango == Procesos - 1){
        for (n = 0; n < N_Species; n++){

            for (j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
                for (k = 0; k < NZ; k++){
                    Species[n].Right[RIGHT(NX,j,k)] = Species[n].Y_Pres[LP(NX-1,j,k,0)];
                }
            }

        }  
    }

    // Here (Null Gradient)
    for (n = 0; n < N_Species; n++){

        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
                Species[n].Here[HERE(i,j,0)] = Species[n].Y_Pres[LP(i,j,0,0)];
            }
        }

    }
    
    // There (Null Gradient)
    for (n = 0; n < N_Species; n++){

        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
                Species[n].There[THERE(i,j,0)] = Species[n].Y_Pres[LP(i,j,NZ-1,0)];
            }
        }

    }

}

// Function to set all the corresponding static temperatures halos (no changing)
void Solver::Get_Species_StaticHalos(){
int i, j, k, n;

    // Bottom (Inlet)
    if (Ix[Rango] <= NX_1){

        for (i = Ix[Rango]; i < min(NX_1, Fx[Rango]); i++){
            for(j = - HP; j < 0; j++){
                for (k = 0; k < NZ; k++){
                    Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Bottom[BOTTOM(i,0,k)];   
                }
            }
        }
        
    }
    
}

// Function to update the temperature halos
void Solver::Get_Species_UpdateHalos(Mesher MESH){
int i, j, k, n;

    // Top (Outlet)
    for (n = 0; n < N_Species; n++){

        for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
            for (j = NY; j < NY + HP; j++){
                for (k = 0; k < NZ; k++){
                     Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Top[TOP(i,NY,k)];
                }
            }
        }

    }

    // Bottom
    if (Fx[Rango] >= NX_1){

        for (i = max(Ix[Rango],NX_1); i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] - HP; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j++){
                for (k = 0; k < NZ; k++){
                    Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Bottom[BOTTOM(i,0,k)];   
                }
            }
        }
        
    }


    // Left
    if (Rango == 0){
        for (n = 0; n < N_Species; n++){

            for (i = - HP; i < 0; i++){
                for (j = 0; j < NY; j++){
                    for (k = 0; k < NZ; k++){
                        Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Left[LEFT(0,j,k)];
                    }
                }
            }
            
        }  
    }

    // Right
    if (Rango == Procesos - 1){
        for (n = 0; n < N_Species; n++){

            for (i = NX; i < NX + HP; i++){
                for (j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
                    for (k = 0; k < NZ; k++){
                        Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Right[RIGHT(NX,j,k)];
                    }
                }
            }
            
        }  
    }

    // Here
    for (n = 0; n < N_Species; n++){

        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
                for(k = - HP; k < 0; k++){
                    Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Here[HERE(i,j,0)];
                }  
            }
        }

    }
    
    // There
    for (n = 0; n < N_Species; n++){

        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
                for (k = NZ; k < NZ + HP; k++){
                    Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].There[THERE(i,j,NZ)];
                }
            }
        }

    }

}