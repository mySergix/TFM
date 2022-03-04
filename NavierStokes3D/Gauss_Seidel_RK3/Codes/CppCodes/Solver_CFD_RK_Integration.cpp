//------------------------------------------------------------------------------------------------//
//                    CPP FILE FOR CFD SOLVER CLASS RUNGE KUTTA INTEGRATION                       //
//------------------------------------------------------------------------------------------------//

// Function to get all the Runge Kutta coefficients
void Solver::Get_RK_Coefficients(Runge_Kutta &RK_Struct){

    // Runge Kutta Coefficients
	RK_Struct.c1 = 0.0;
	RK_Struct.c2 = 0.50;
	RK_Struct.c3 = 1.0;

	RK_Struct.a11 = 0.0;
	RK_Struct.a12 = 0.0;
	RK_Struct.a13 = 0.0;

	RK_Struct.a21 = 0.50;
	RK_Struct.a22 = 0.0;
	RK_Struct.a23 = 0.0;

	RK_Struct.a31 = - 1.0;
	RK_Struct.a32 = 2.0;
	RK_Struct.a33 = 0.0;

	RK_Struct.b1 = 1.0 / 6.0;
	RK_Struct.b2 = 2.0 / 3.0;
	RK_Struct.b3 = 1.0 / 6.0;

}

// Function to sum the contributions of the flow
void Solver::Get_RK_VelocityContributions(double *Contribution_Matrix_U, double *Contribution_Matrix_V, double *Contribution_Matrix_W){
int i, j, k;

    // Velocity U
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				Contribution_Matrix_U[LU(i,j,k,0)] = U.Diffusive[LU(i,j,k,0)] - U.Convective[LU(i,j,k,0)];
			}
		}
	}

	// Velocity V
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 1; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				Contribution_Matrix_V[LV(i,j,k,0)] = V.Diffusive[LV(i,j,k,0)] - V.Convective[LV(i,j,k,0)] + V.Boussinesq[LV(i,j,k,0)];
			}
		}
	}

    // Velocity W
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 1; k < NZ; k++){
				Contribution_Matrix_W[LW(i,j,k,0)] = W.Diffusive[LW(i,j,k,0)] - W.Convective[LW(i,j,k,0)];
			}
		}
	}

}

// Function to get the Runge Kutta temporal integration
void Solver::Get_RK3_Integration(Mesher MESH, Parallel P1){
int i, j, k;

    // 1. Calculo Diffusion + Convection + Boussinesq
        
        // Boundary Conditions Update
        Get_UpdateBoundaryConditions_Velocities(U.Pres, V.Pres, W.Pres);
        Get_UpdateHalos_Velocity(U.Pres, V.Pres, W.Pres);

        // Halos communication
		P1.CommunicateDataLU(U.Pres, U.Pres);
		P1.CommunicateDataLV(V.Pres, V.Pres);
		P1.CommunicateDataLW(W.Pres, W.Pres);
        
        // Diffusion Terms Calculation
        Get_DiffusionU(MESH, U.Pres);
        Get_DiffusionV(MESH, V.Pres);
        Get_DiffusionW(MESH, W.Pres);
        
		// Convective Terms Calculation
		Get_ConvectionU(MESH, U.Pres, V.Pres, W.Pres);
		Get_ConvectionV(MESH, U.Pres, V.Pres, W.Pres);
		Get_ConvectionW(MESH, U.Pres, V.Pres, W.Pres);

        // Energy Equation
        if (Problema == 2){	
			Get_UpdateBoundaryConditions_Temperatures(T.Pres);
			Get_UpdateHalos_Temperatures(T.Pres);

            P1.CommunicateDataLP(T.Pres, T.Pres); 
            Get_DiffusionEnergy(MESH, T.Pres);
			Get_ConvectionEnergy(MESH, T.Pres, U.Pres, V.Pres, W.Pres);
            Get_RK_ScalarContributions(T.K1, T.Diffusive, T.Convective);
            Get_BoussinesqV(MESH, V.Pres, T.Pres); 
        }

    // 2. Calculo K1    
    
        Get_RK_VelocityContributions(U.K1, V.K1, W.K1);
        
    // 3. Calculo Nueva Velocidad

        // Intermediate Velocity Predictor U
        for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
            for(j = 0; j < NY; j++){
		        for(k = 0; k < NZ; k++){
				    U.New_Velocity[LU(i,j,k,0)] = U.Pres[LU(i,j,k,0)] + RK3.a21 * DeltaT * U.K1[LU(i,j,k,0)];
			    }
		    }
	    }

        // Intermediate Velocity Predictor V
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
		    for(k = 0; k < NZ; k++){

			    // Bottom Part
			    V.New_Velocity[LV(i,0,k,0)] = V.Bottom[BOTTOM(i,0,k)];

        	    for(j = 1; j < NY; j++){
				    V.New_Velocity[LV(i,j,k,0)] = V.Pres[LV(i,j,k,0)] + RK3.a21 * DeltaT * V.K1[LV(i,j,k,0)];
			    }

			    // Top Part
			    V.New_Velocity[LV(i,NY,k,0)] = V.Top[TOP(i,NY,k)];

		    }
	    }

        // Intermediate Velocity Predictor W
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = 0; j < NY; j++){

			    // Here Part
			    W.New_Velocity[LW(i,j,0,0)] = W.Here[HERE(i,j,0)];
			
		        for(k = 1; k < NZ; k++){
				    W.New_Velocity[LW(i,j,k,0)] = W.Pres[LW(i,j,k,0)] + RK3.a21 * DeltaT * W.K1[LW(i,j,k,0)];
			    }

			    // There Part
			    W.New_Velocity[LW(i,j,0,0)] = W.There[THERE(i,j,NZ)];

		    }
	    }

        // Calculo nueva temperatura
        if (Problema == 2){

            // Intermediate Temperature Calculation
            for(i = Ix[Rango]; i < Fx[Rango]; i++){
                for(j = 0; j < NY; j++){
		            for(k = 0; k < NZ; k++){
				        T.New_Temperature[LP(i,j,k,0)] = T.Pres[LP(i,j,k,0)] + RK3.a21 * DeltaT * T.K1[LP(i,j,k,0)];
			        }
		        }
	        }

        }

    // 4. Calculo Diffusion + Convection + Boussinesq

        // Boundary Conditions Update
        Get_UpdateBoundaryConditions_Velocities(U.New_Velocity, V.New_Velocity, W.New_Velocity);
        Get_UpdateHalos_Velocity(U.New_Velocity, V.New_Velocity, W.New_Velocity);

        // Halos communication
		P1.CommunicateDataLU(U.New_Velocity, U.New_Velocity);
		P1.CommunicateDataLV(V.New_Velocity, V.New_Velocity);
		P1.CommunicateDataLW(W.New_Velocity, W.New_Velocity);

        // Diffusion Terms Calculation
        Get_DiffusionU(MESH, U.New_Velocity);
        Get_DiffusionV(MESH, V.New_Velocity);
        Get_DiffusionW(MESH, W.New_Velocity);
        
		// Convective Terms Calculation
		Get_ConvectionU(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);
		Get_ConvectionV(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);
		Get_ConvectionW(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);

        // Energy Equation
        if (Problema == 2){	
			Get_UpdateBoundaryConditions_Temperatures(T.New_Temperature);
			Get_UpdateHalos_Temperatures(T.New_Temperature);

            P1.CommunicateDataLP(T.New_Temperature, T.New_Temperature); 
            Get_DiffusionEnergy(MESH, T.New_Temperature);
			Get_ConvectionEnergy(MESH, T.New_Temperature, U.New_Velocity, V.New_Velocity, W.New_Velocity);
            Get_RK_ScalarContributions(T.K2, T.Diffusive, T.Convective);
            Get_BoussinesqV(MESH, V.New_Velocity, T.New_Temperature); 
        }

    // 5. Calculo K2;

        Get_RK_VelocityContributions(U.K2, V.K2, W.K2);

    // 6. Calculo Nueva Velocidad

        // Intermediate Velocity Predictor U
        for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
            for(j = 0; j < NY; j++){
		        for(k = 0; k < NZ; k++){
				    U.New_Velocity[LU(i,j,k,0)] = U.Pres[LU(i,j,k,0)] + DeltaT * (RK3.a31 * U.K1[LU(i,j,k,0)] + RK3.a32 * U.K2[LU(i,j,k,0)]);
			    }
		    }
	    }

        // Intermediate Velocity Predictor V
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
		    for(k = 0; k < NZ; k++){

			    // Bottom Part
			    V.New_Velocity[LV(i,0,k,0)] = V.Bottom[BOTTOM(i,0,k)];

        	    for(j = 1; j < NY; j++){
				    V.New_Velocity[LV(i,j,k,0)] = V.Pres[LV(i,j,k,0)] + DeltaT * (RK3.a31 * V.K1[LV(i,j,k,0)] + RK3.a32 * V.K2[LV(i,j,k,0)]);
			    }

			    // Top Part
			    V.New_Velocity[LV(i,NY,k,0)] = V.Top[TOP(i,NY,k)];

		    }
	    }

        // Intermediate Velocity Predictor W
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = 0; j < NY; j++){

			    // Here Part
			    W.New_Velocity[LW(i,j,0,0)] = W.Here[HERE(i,j,0)];
			
		        for(k = 1; k < NZ; k++){
				    W.New_Velocity[LW(i,j,k,0)] = W.Pres[LW(i,j,k,0)] + DeltaT * (RK3.a31 * W.K1[LW(i,j,k,0)] + RK3.a32 * W.K2[LW(i,j,k,0)]);
			    }

			    // There Part
			    W.New_Velocity[LW(i,j,0,0)] = W.There[THERE(i,j,NZ)];

		    }
	    }

        // Calculo nueva temperatura
        if (Problema == 2){

            // Intermediate Temperature Calculation
            for(i = Ix[Rango]; i < Fx[Rango]; i++){
                for(j = 0; j < NY; j++){
		            for(k = 0; k < NZ; k++){
				        T.New_Temperature[LP(i,j,k,0)] = T.Pres[LP(i,j,k,0)] + DeltaT * (RK3.a31 * T.K1[LP(i,j,k,0)] + RK3.a32 * T.K2[LP(i,j,k,0)]);
			        }
		        }
	        }

        }

    // 7. Calculo Diffusion + Convection + Boussinesq

        // Boundary Conditions Update
        Get_UpdateBoundaryConditions_Velocities(U.New_Velocity, V.New_Velocity, W.New_Velocity);
        Get_UpdateHalos_Velocity(U.New_Velocity, V.New_Velocity, W.New_Velocity);

        // Halos communication
		P1.CommunicateDataLU(U.New_Velocity, U.New_Velocity);
		P1.CommunicateDataLV(V.New_Velocity, V.New_Velocity);
		P1.CommunicateDataLW(W.New_Velocity, W.New_Velocity);

        // Diffusion Terms Calculation
        Get_DiffusionU(MESH, U.New_Velocity);
        Get_DiffusionV(MESH, V.New_Velocity);
        Get_DiffusionW(MESH, W.New_Velocity);
        
		// Convective Terms Calculation
		Get_ConvectionU(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);
		Get_ConvectionV(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);
		Get_ConvectionW(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);

        // Energy Equation
        if (Problema == 2){	
			Get_UpdateBoundaryConditions_Temperatures(T.New_Temperature);
			Get_UpdateHalos_Temperatures(T.New_Temperature);

            P1.CommunicateDataLP(T.New_Temperature, T.New_Temperature); 
            Get_DiffusionEnergy(MESH, T.New_Temperature);
			Get_ConvectionEnergy(MESH, T.New_Temperature, U.New_Velocity, V.New_Velocity, W.New_Velocity);
            Get_RK_ScalarContributions(T.K3, T.Diffusive, T.Convective);
            Get_BoussinesqV(MESH, V.New_Velocity, T.New_Temperature); 
        }

    // 8. Calculo K3;

        Get_RK_VelocityContributions(U.K3, V.K3, W.K3);

    // 9. Calculo predictora final

        // Velocity Predictor U
        for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
            for(j = 0; j < NY; j++){
		        for(k = 0; k < NZ; k++){
				    U.Predictor[LU(i,j,k,0)] = U.Pres[LU(i,j,k,0)] + DeltaT * (RK3.b1 * U.K1[LU(i,j,k,0)] + RK3.b2 * U.K2[LU(i,j,k,0)] + RK3.b3 * U.K3[LU(i,j,k,0)]);
			    }
		    }
	    }

        // HAY QUE ACTUALIZAR LAS CONDICIONES DE CONTORNO AQUI
        // Velocity Predictor V
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
		    for(k = 0; k < NZ; k++){

			    // Bottom Part
			    V.Predictor[LV(i,0,k,0)] = V.Bottom[BOTTOM(i,0,k)];

        	    for(j = 1; j < NY; j++){
				    V.Predictor[LV(i,j,k,0)] = V.Pres[LV(i,j,k,0)] + DeltaT * (RK3.b1 * V.K1[LV(i,j,k,0)] + RK3.b2 * V.K2[LV(i,j,k,0)] + RK3.b3 * V.K3[LV(i,j,k,0)]);
			    }

			    // Top Part
			    V.Predictor[LV(i,NY,k,0)] = V.Top[TOP(i,NY,k)];

		    }
	    }

        // Velocity Predictor W
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = 0; j < NY; j++){

			    // Here Part
			    W.Predictor[LW(i,j,0,0)] = W.Here[HERE(i,j,0)];
			
		        for(k = 1; k < NZ; k++){
				    W.Predictor[LW(i,j,k,0)] = W.Pres[LW(i,j,k,0)] + DeltaT * (RK3.b1 * W.K1[LW(i,j,k,0)] + RK3.b2 * W.K2[LW(i,j,k,0)] + RK3.b3 * W.K3[LW(i,j,k,0)]);
			    }

			    // There Part
			    W.Predictor[LW(i,j,NZ,0)] = W.There[THERE(i,j,NZ)];

		    }
	    }

        Get_UpdateBoundaryConditions_PredictorVelocities();

        // Calculo nueva temperatura
        if (Problema == 2){

            // Final Temperature Calculation
            for(i = Ix[Rango]; i < Fx[Rango]; i++){
                for(j = 0; j < NY; j++){
		            for(k = 0; k < NZ; k++){
				        T.Fut[LP(i,j,k,0)] = T.Pres[LP(i,j,k,0)] + DeltaT * (RK3.b1 * T.K1[LP(i,j,k,0)] + RK3.b2 * T.K2[LP(i,j,k,0)] + RK3.b3 * T.K3[LP(i,j,k,0)]);
			        }
		        }
	        }

        }

}

// Function to sum the contributions of the scalar property
void Solver::Get_RK_ScalarContributions(double *Contribution_Matrix, double *Contribution_Diffusive, double *Contribution_Convective){
int i, j, k;

	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				Contribution_Matrix[LP(i,j,k,0)] = Contribution_Diffusive[LP(i,j,k,0)] - Contribution_Convective[LP(i,j,k,0)];
			}
		}
	}

}