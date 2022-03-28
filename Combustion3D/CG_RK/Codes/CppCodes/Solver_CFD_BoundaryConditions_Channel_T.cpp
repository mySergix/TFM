//------------------------------------------------------------------------------------------------//
//               CPP FILE FOR MOIST CHANNEL TEMPERATURE BOUNDARY CONDITIONS SETUP                 //
//------------------------------------------------------------------------------------------------//

// Function to set all the initial Velocity Boundary Conditions of the simulation
void Solver::Get_StaticBoundaryConditions_Temperatures(Mesher MESH){
int i, j, k;

    // Left (Inlet)
    if (Rango == 0){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                T.Left[LEFT(0,j,k)] = To;
            }
        }
    }

    // Bottom
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            if (MESH.MP[LP(i,0,k,2)] >= 0.009 && MESH.MP[LP(i,0,k,2)] <= 0.289){
                T.Bottom[BOTTOM(i,0,k)] = Twater;
            }
        }
    }
    
}

// Function to update the temperatures boundary conditions
void Solver::Get_UpdateBoundaryConditions_Temperatures(Mesher MESH, double *T_Matrix){
int i, j, k, n;

    // Species 0 -> Water Vapour
    n = 0;
    
    // Right (Outlet) (Null Gradient)
    if (Rango == Procesos - 1){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                T.Right[RIGHT(NX,j,k)] = T_Matrix[LP(NX-1,j,k,0)];
            }
        }
    }  

    // Here
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            T.Here[HERE(i,j,0)] = T_Matrix[LP(i,j,0,0)];
        }
    }

    // There
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            T.There[THERE(i,j,NZ)] = T_Matrix[LP(i,j,NZ-1,0)];
        }
    }

    // Top
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
            T.Top[TOP(i,NY,k)] = T_Matrix[LP(i,NY-1,k,0)];
        }
    }

    // Bottom
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            if (MESH.MP[LP(i,0,k,2)] >= 0.009 && MESH.MP[LP(i,0,k,2)] <= 0.289){
                T.Bottom[BOTTOM(i,0,k)]  = T_Matrix[LP(i,0,k,0)] + (D_AB / K) * (Species[n].Y_Pres[LP(i,0,k,0)] - Species[n].Bottom[BOTTOM(i,0,k)]) * hfg;
                
            }
            else{
                T.Bottom[BOTTOM(i,0,k)] = T_Matrix[LP(i,0,k,0)];
            }
        }
    }

}

// Function to set all the corresponding static temperatures halos (no changing)
void Solver::Get_StaticHalos_Temperatures(double *T_Matrix){
int i, j, k;

    // Left (Inlet)
    if (Rango == 0){
        for (i = - Halo; i < Ix[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    T_Matrix[LP(i,j,k,0)] = T.Left[LEFT(0,j,k)];   
                }
            }
        }  
    }
    
    // Bottom
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = - Halo; j < 0; j++){
            for (k = 0; k < NZ; k++){
                T_Matrix[LP(i,j,k,0)] = T.Bottom[BOTTOM(i,0,k)];
            }
        }   
    }

}

// Function to update the temperature halos
void Solver::Get_UpdateHalos_Temperatures(double *T_Matrix){
int i, j, k;

    // Right (Outlet)
    if (Rango == Procesos - 1){
        for (i = NX; i < NX + Halo; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    T_Matrix[LP(i,j,k,0)] = T.Right[RIGHT(NX,j,k)];
                }
            }
        }    
    }  

    // Top
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = NY; j < NY + Halo; j++){
            for (k = 0; k < NZ; k++){
                T_Matrix[LP(i,j,k,0)] = T.Top[TOP(i,NY,k)];
            }
        }  
    }

    // Here
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for(k = - Halo; k < 0; k++){
                T_Matrix[LP(i,j,k,0)] = T.Here[HERE(i,j,0)];
            }  
        }
    }

    // There
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = NZ; k < NZ + Halo; k++){
                T_Matrix[LP(i,j,k,0)] = T.There[THERE(i,j,NZ)];
            }
        }
    }

}