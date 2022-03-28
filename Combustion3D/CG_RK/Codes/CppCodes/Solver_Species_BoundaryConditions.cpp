//------------------------------------------------------------------------------------------------//
//                        CPP FILE FOR SPECIES BOUNDARY CONDITIONS SETUP                          //
//------------------------------------------------------------------------------------------------//

// Function to set all the initial Species Boundary Conditions of the simulation
void Solver::Get_Species_StaticBoundaryConditions(Mesher MESH){
int i, j, k, n;

    // Species 0 -> Water Vapour
    n = 0;

    // Left (Inlet)
    if (Rango == 0){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[n].Left[LEFT(0,j,k)] = Co;
            }
        }
    }

}

// Function to update the temperatures boundary conditions
void Solver::Get_Species_UpdateBoundaryConditions(Mesher MESH){
int i, j, k, n;

    // Species 0 -> Water Vapour
    n = 0;

    // Right (Outlet) (Null Gradient)
    if (Rango == Procesos - 1){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[n].Right[RIGHT(NX,j,k)] = Species[n].Y_Pres[LP(NX-1,j,k,0)];
            }
        }
    }  

    // Here (Null Gradient)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            Species[n].Here[HERE(i,j,0)] = Species[n].Y_Pres[LP(i,j,0,0)];
        }
    }

    // There (Null Gradient)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            Species[n].There[THERE(i,j,NZ)] = Species[n].Y_Pres[LP(i,j,NZ-1,0)];
        }
    }

    // Bottom (Fixed Value)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
            if (MESH.MP[LP(i,0,k,2)] >= 0.009 && MESH.MP[LP(i,0,k,2)] <= 0.289){ 
                Species[n].Bottom[BOTTOM(i,0,k)] = (100 * pow(10, 0.66077 + (7.5 * T.Bottom[BOTTOM(i,0,k)]) / (237.3 + T.Bottom[BOTTOM(i,0,k)])) * (MW_Water / 1000)) / (8.314 * (273.15 + T.Bottom[BOTTOM(i,0,k)]));
            }
            else{
                Species[n].Bottom[BOTTOM(i,0,k)] = Species[n].Y_Pres[LP(i,0,k,0)];
            }
        }
    }

    // Top (Null Gradient)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
            Species[n].Top[TOP(i,NY,k)] = Species[n].Y_Pres[LP(i,NY-1,k,0)];
        }
    }

}

// Function to set all the corresponding static temperatures halos (no changing)
void Solver::Get_Species_StaticHalos(){
int i, j, k, n;

    // Species 0 -> Water Vapour
    n = 0;

    // Left (Inlet)
    if (Rango == 0){
        for (i = - Halo; i < Ix[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Left[LEFT(0,j,k)];   
                }
            }
        }  
    }

}

// Function to update the temperature halos
void Solver::Get_Species_UpdateHalos(){
int i, j, k, n;

    // Species 0 -> Water Vapour
    n = 0;

    // Right (Outlet)
    if (Rango == Procesos - 1){
        for (i = NX; i < NX + Halo; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Right[RIGHT(NX,j,k)];
                }
            }
        }    
    }  

    // Bottom
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = - Halo; j < 0; j++){
            for (k = 0; k < NZ; k++){
                Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Bottom[BOTTOM(i,0,k)];
            }
        }  
    }

    // Top
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = NY; j < NY + Halo; j++){
            for (k = 0; k < NZ; k++){
                Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Top[TOP(i,NY,k)];
            }
        }  
    }

    // Here
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for(k = - Halo; k < 0; k++){
                Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Here[HERE(i,j,0)];
            }  
        }
    }

    // There
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = NZ; k < NZ + Halo; k++){
                Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].There[THERE(i,j,NZ)];
            }
        }
    }

}