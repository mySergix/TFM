//------------------------------------------------------------------------------------------------//
//               CPP FILE FOR MOIST CHANNEL VELOCITY BOUNDARY CONDITIONS SETUP                    //
//------------------------------------------------------------------------------------------------//

// Function to set all the initial Velocity Boundary Conditions of the simulation
void Solver::Get_StaticBoundaryConditions_Velocities(Mesher MESH){
int i, j, k;

    // Left (Inlet)
    double m = 1.7 + 0.50 * pow(Gamma_Geometry, -1.4);
    double n = 2;// + 0.3 * (Gamma_Geometry - 1/3);

    if (Rango == 0){ // Left
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Left[LEFT(0,j,k)] = w_av * ((m+1)/m) * ((n+1)/n) * (1.0 - pow(abs(MESH.MP[LP(0,j,k,1)] - 0.50 * Ydominio) / (0.50 * Ydominio), n)) * (1.0 - pow(abs(MESH.MP[LP(0,j,k,2)] - 0.50 * Zdominio) / (0.50 * Zdominio), m));
                V.Left[LEFT(0,j,k)] = 0.0;
                W.Left[LEFT(0,j,k)] = 0.0;
            }
        }
    }
    
    // Bottom
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            U.Bottom[BOTTOM(i,0,k)] = 0.0;
            V.Bottom[BOTTOM(i,0,k)] = 0.0;
            W.Bottom[BOTTOM(i,0,k)] = 0.0;
        }
    }

    // Top
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (k = 0; k < NZ; k++){
            U.Top[TOP(i,NY,k)] = 0.0;
            V.Top[TOP(i,NY,k)] = 0.0;
            W.Top[TOP(i,NY,k)] = 0.0;
        }
    }
    
    // Here
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            U.Here[HERE(i,j,0)] = 0.0;
            V.Here[HERE(i,j,0)] = 0.0;
            W.Here[HERE(i,j,0)] = 0.0;
        }
    }

    // There
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            U.There[THERE(i,j,NZ)] = 0.0;
            V.There[THERE(i,j,NZ)] = 0.0;
            W.There[THERE(i,j,NZ)] = 0.0;
        }
    }

}

// Function to update the velocity boundary conditions
void Solver::Get_UpdateBoundaryConditions_Velocities(double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;

    // Right (Outlet) (Null Gradient)
    if (Rango == Procesos - 1){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Right[RIGHT(NX,j,k)] = U_Matrix[LU(NX-1,j,k,0)];
                V.Right[RIGHT(NX,j,k)] = V_Matrix[LV(NX-1,j,k,0)];
                W.Right[RIGHT(NX,j,k)] = W_Matrix[LW(NX-1,j,k,0)];
            }
        }
    }
    
}

// Function to update the predictor velocities boundary conditions
void Solver::Get_UpdateBoundaryConditions_PredictorVelocities(){
int i, j, k;

    // Left (Inlet)
    if (Rango == 0){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Predictor[LU(0,j,k,0)] = U.Left[LEFT(0,j,k)];
            }
        }
    }

    // Right (Outlet)
    if (Rango == Procesos - 1){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Predictor[LU(NX,j,k,0)] = U.Predictor[LU(NX-1,j,k,0)];
            }
        }
    }

    // Bottom (Non-Slip)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
			V.Predictor[LV(i,0,k,0)] = V.Pres[LV(i,0,k,0)];
        }
    }

    // Top (Non-Slip)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
			V.Predictor[LV(i,NY,k,0)] = V.Pres[LV(i,NY,k,0)];
        }
    }

    // Here (Non-Slip)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            W.Predictor[LW(i,j,0,0)] = W.Here[HERE(i,j,0)];
        }
    }

    // There (Non-Slip)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            W.Predictor[LW(i,j,NZ,0)] = W.There[THERE(i,j,NZ)];
        }
    }

}

// Function to set all the corresponding static velocity halos (no changing)
void Solver::Get_StaticHalos_Velocity(double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;

    // Left (Inlet)
    if (Rango == 0){ 
        for (i = - Halo; i < Ix[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    V_Matrix[LV(i,j,k,0)] = V.Left[LEFT(0,j,k)];
                    W_Matrix[LW(i,j,k,0)] = W.Left[LEFT(0,j,k)];
                }
            }
        }  

        for (i = - Halo; i <= Ix[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    U_Matrix[LU(i,j,k,0)] = U.Left[LEFT(0,j,k)];
                }
            }
        }  
    }
    
    // Bottom
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = - Halo; j < 0; j++){
            for (k = 0; k < NZ; k++){
                U_Matrix[LU(i,j,k,0)] = U.Bottom[BOTTOM(i,0,k)];
                W_Matrix[LW(i,j,k,0)] = W.Bottom[BOTTOM(i,0,k)];
            }
        }   
    }
    
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = - Halo; j <= 0; j++){
            for (k = 0; k < NZ; k++){
                V_Matrix[LV(i,j,k,0)] = V.Bottom[BOTTOM(i,0,k)];
            }
        }   
    }

    // Top
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = NY; j < NY + Halo; j++){
            for (k = 0; k < NZ; k++){
                U_Matrix[LU(i,j,k,0)] = U.Top[TOP(i,NY,k)];
                W_Matrix[LW(i,j,k,0)] = W.Top[TOP(i,NY,k)];
            }
        }  
    }

    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = NY; j < NY + Halo + 1; j++){
            for (k = 0; k < NZ; k++){
                V_Matrix[LV(i,j,k,0)] = V.Top[TOP(i,NY,k)];
            }
        }  
    }
 
    // Here
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for(k = - Halo; k < 0; k++){
                U_Matrix[LU(i,j,k,0)] = U.Here[HERE(i,j,0)];
                V_Matrix[LV(i,j,k,0)] = V.Here[HERE(i,j,0)];
            }  
        }
    }

    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for(k = - Halo; k <= 0; k++){
                W_Matrix[LW(i,j,k,0)] = W.Here[HERE(i,j,0)];
            }  
        }
    }


    // There
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = NZ; k < NZ + Halo; k++){
                U_Matrix[LU(i,j,k,0)] = U.There[THERE(i,j,NZ)];
                V_Matrix[LV(i,j,k,0)] = V.There[THERE(i,j,NZ)];
            }    
        }
    }

    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = NZ; k < NZ + Halo + 1; k++){
                W_Matrix[LW(i,j,k,0)] = W.There[THERE(i,j,NZ)];
            }     
        }
    }

}

// Function to update the boundary conditions
void Solver::Get_UpdateHalos_Velocity(double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;

    // Right (Outlet)
    if (Rango == Procesos - 1){
        for (i = Fx[Rango]; i < Fx[Rango] + Halo; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    V_Matrix[LV(i,j,k,0)] = V.Right[RIGHT(NX,j,k)];
                    W_Matrix[LW(i,j,k,0)] = W.Right[RIGHT(NX,j,k)];
                }
            }
        }   

        for (i = Fx[Rango]; i < Fx[Rango] + Halo + 1; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    U_Matrix[LU(i,j,k,0)] = U.Right[RIGHT(NX,j,k)];
                }
            }
        } 
    } 

}