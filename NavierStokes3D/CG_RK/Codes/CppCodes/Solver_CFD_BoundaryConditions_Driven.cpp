//------------------------------------------------------------------------------------------------//
//                      CPP FILE FOR DRIVEN CAVITY BOUNDARY CONDITIONS SETUP                      //
//------------------------------------------------------------------------------------------------//

// Function to set all the initial Boundary Conditions of the simulation
void Solver::Get_StaticBoundaryConditions_Velocities(Mesher MESH){
int i, j, k;

    // Left
    if (Rango == 0){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Left[LEFT(0,j,k)] = 0.0;
                V.Left[LEFT(0,j,k)] = 0.0;
                W.Left[LEFT(0,j,k)] = 0.0;
            }
        }
    }

    // Right
    if (Rango == Procesos - 1){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Right[RIGHT(NX,j,k)] = 0.0;
                V.Right[RIGHT(NX,j,k)] = 0.0;
                W.Right[RIGHT(NX,j,k)] = 0.0;
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
            U.Top[TOP(i,NY,k)] = Uref;
            V.Top[TOP(i,NY,k)] = 0.0;
            W.Top[TOP(i,NY,k)] = 0.0;
        }
    }
    
}

// Function to update the velocity boundary conditions
void Solver::Get_UpdateBoundaryConditions_Velocities(double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;

    // Here
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            U.Here[HERE(i,j,0)] = U_Matrix[LU(i,j,0,0)];
            V.Here[HERE(i,j,0)] = V_Matrix[LV(i,j,0,0)];
        }
    }

    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            W.Here[HERE(i,j,0)] = W_Matrix[LW(i,j,1,0)];
        }
    }


    // There
    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            U.There[THERE(i,j,NZ)] = U_Matrix[LU(i,j,NZ-1,0)];
            V.There[THERE(i,j,NZ)] = V_Matrix[LV(i,j,NZ-1,0)];
        }
    }

    for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            W.There[THERE(i,j,NZ)] = W_Matrix[LW(i,j,NZ-1,0)];
        }
    }

}

// Function to update the predictor velocities boundary conditions
void Solver::Get_UpdateBoundaryConditions_PredictorVelocities(){
int i, j, k;

    // Left (Non-Slip)
    if (Rango == 0){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Predictor[LU(0,j,k,0)] = U.Pres[LU(0,j,k,0)];
            }
        }
    }

    // Right (Non-Slip)
    if (Rango == Procesos - 1){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                U.Predictor[LU(NX,j,k,0)] = U.Pres[LU(NX,j,k,0)];
            }
        }
    }

    // Bottom (Non-Slip)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
			V.Predictor[LV(i,0,k,0)] = V.Pres[LV(i,0,k,0)];
        }
    }

    // Top (V = 0)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 0; k < NZ; k++){
			V.Predictor[LV(i,NY,k,0)] = V.Pres[LV(i,NY,k,0)];
        }
    }

    // Here (Null Gradient)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            W.Predictor[LW(i,j,0,0)] = W.Predictor[LW(i,j,1,0)];
        }
    }

    // There (Null Gradient)
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            W.Predictor[LW(i,j,NZ,0)] = W.Predictor[LW(i,j,NZ-1,0)];
        }
    }

}

// Function to set all the corresponding static velocity halos (no changing)
void Solver::Get_StaticHalos_Velocity(double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;

    // Left
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

    // Right
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
 
}

// Function to update the boundary conditions
void Solver::Get_UpdateHalos_Velocity(double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;

    // Here
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for(k = - Halo; k < 0; k++){
                U_Matrix[LU(i,j,k,0)] = U.Here[HERE(i,j,0)];
                V_Matrix[LV(i,j,k,0)] = V.Here[HERE(i,j,0)];
            }  
        }
    }

    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for(k = - Halo; k <= 0; k++){
                W_Matrix[LW(i,j,k,0)] = W.Here[HERE(i,j,0)];
            }  
        }
    }

    // There
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = NZ; k < NZ + Halo; k++){
                U_Matrix[LU(i,j,k,0)] = U.There[THERE(i,j,NZ)];
                V_Matrix[LV(i,j,k,0)] = V.There[THERE(i,j,NZ)];
            } 
        }
    }

    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (j = 0; j < NY; j++){
            for (k = NZ; k < NZ + Halo + 1; k++){
                W_Matrix[LW(i,j,k,0)] = W.There[THERE(i,j,NZ)];
            } 
        }
    }

}