//------------------------------------------------------------------------------------------------//
//                    CPP FILE FOR CFD SOLVER CLASS RUNGE KUTTA 4 INTEGRATION                     //
//------------------------------------------------------------------------------------------------//

// Function to get all the Runge Kutta coefficients
void Solver::Get_RK_Coefficients(Runge_Kutta &RK_Struct){

    // Runge Kutta Coefficients
	RK_Struct.c1 = 0.0;
	RK_Struct.c2 = 0.50;
	RK_Struct.c3 = 0.50;
	RK_Struct.c4 = 1.0;

	RK_Struct.a11 = 0.0;
	RK_Struct.a12 = 0.0;
	RK_Struct.a13 = 0.0;
	RK_Struct.a14 = 0.0;

	RK_Struct.a21 = 0.50;
	RK_Struct.a22 = 0.0;
	RK_Struct.a23 = 0.0;
	RK_Struct.a24 = 0.0;

	RK_Struct.a31 = 0.0;
	RK_Struct.a32 = 0.50;
	RK_Struct.a33 = 0.0;
	RK_Struct.a34 = 0.0;

	RK_Struct.a41 = 0.0;
	RK_Struct.a42 = 0.0;
	RK_Struct.a43 = 1.0;
	RK_Struct.a44 = 0.0;

	RK_Struct.b1 = 1.0 / 6.0;
	RK_Struct.b2 = 2.0 / 6.0;
	RK_Struct.b3 = 2.0 / 6.0;
	RK_Struct.b4 = 1.0 / 6.0;

}

// Function to sum the contributions of the flow
void Solver::Get_RK_VelocityContributions(Mesher MESH, double *Contribution_Matrix_U, double *Contribution_Matrix_V, double *Contribution_Matrix_W){
int i, j, k;

    // Velocity U
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1]; j++){
		    for(k = 0; k < NZ; k++){
				Contribution_Matrix_U[LU(i,j,k,0)] = DeltaT * (U.Diffusive[LU(i,j,k,0)] - U.Convective[LU(i,j,k,0)]);
			}
		}
	}

	// Velocity V
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){
		    for(k = 0; k < NZ; k++){
				Contribution_Matrix_V[LV(i,j,k,0)] = DeltaT * (V.Diffusive[LV(i,j,k,0)] - V.Convective[LV(i,j,k,0)] + V.Boussinesq[LV(i,j,k,0)]);
			}
		}
	}

    // Velocity W
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMW[i + Halo - Ix[Rango]][1]; j++){
		    for(k = 1; k < NZ; k++){
				Contribution_Matrix_W[LW(i,j,k,0)] = DeltaT * (W.Diffusive[LW(i,j,k,0)] - W.Convective[LW(i,j,k,0)]);
			}
		}
	}

}

// Function to get the Runge Kutta temporal integration
void Solver::Get_RK_Integration(Mesher MESH, Parallel P1){
int i, j, k;

    // 1. Calculo Diffusion + Convection + Boussinesq
        
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

        // Boussinesq Thermal Buoyancy Term
        //Get_BoussinesqV(MESH, V.Pres, T.Pres);

    // 2. Calculo K1    
    
        Get_RK_VelocityContributions(MESH, U.K1, V.K1, W.K1);
    
    // 3. Calculo Nueva Velocidad
	
        // Intermediate Velocity Predictor U
        for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
            for(j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1]; j++){
		        for(k = 0; k < NZ; k++){
				    U.New_Velocity[LU(i,j,k,0)] = U.Pres[LU(i,j,k,0)] + RK.a21 * U.K1[LU(i,j,k,0)];
			    }
		    }
	    }

        // Intermediate Velocity Predictor V
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
			for(j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){
		    	for(k = 0; k < NZ; k++){  
				    V.New_Velocity[LV(i,j,k,0)] = V.Pres[LV(i,j,k,0)] + RK.a21 * V.K1[LV(i,j,k,0)];
			    }
		    }
	    }

        // Intermediate Velocity Predictor W
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
			for(j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMW[i + Halo - Ix[Rango]][1]; j++){
		        for(k = 1; k < NZ; k++){
				    W.New_Velocity[LW(i,j,k,0)] = W.Pres[LW(i,j,k,0)] + RK.a21 * W.K1[LW(i,j,k,0)];
			    }
		    }
	    }

	// 4. Intermediate Velocities Correction
	
		Get_CorrectedVelocities(MESH, P1, U.New_Velocity, V.New_Velocity, W.New_Velocity);

    // 5. Calculo Diffusion + Convection + Boussinesq

        // Boundary Conditions Update
        Get_UpdateBoundaryConditions_Velocities(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);
        Get_UpdateHalos_Velocity(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);

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

        // Boussinesq Thermal Buoyancy Term
        //Get_BoussinesqV(MESH, V.New_Velocity, T.Pres);

    // 6. Calculo K2;

        Get_RK_VelocityContributions(MESH, U.K2, V.K2, W.K2);
	
    // 7. Calculo Nueva Velocidad

        // Intermediate Velocity Predictor U
        for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
            for(j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1]; j++){
		        for(k = 0; k < NZ; k++){
				    U.New_Velocity[LU(i,j,k,0)] = U.Pres[LU(i,j,k,0)] + RK.a32 * U.K2[LU(i,j,k,0)];
			    }
		    }
	    }

        // Intermediate Velocity Predictor V
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
			for(j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1] - 1; j++){
		    	for(k = 0; k < NZ; k++){ 
				    V.New_Velocity[LV(i,j,k,0)] = V.Pres[LV(i,j,k,0)] + RK.a32 * V.K2[LV(i,j,k,0)];
			    }
		    }
	    }

        // Intermediate Velocity Predictor W
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMW[i + Halo - Ix[Rango]][1]; j++){
		        for(k = 1; k < NZ; k++){
				    W.New_Velocity[LW(i,j,k,0)] = W.Pres[LW(i,j,k,0)] + RK.a32 * W.K2[LW(i,j,k,0)];
			    }
		    }
	    }

	// 8. Intermediate Velocities Correction
	
		Get_CorrectedVelocities(MESH, P1, U.New_Velocity, V.New_Velocity, W.New_Velocity);

    // 9. Calculo Diffusion + Convection + Boussinesq

        // Boundary Conditions Update
        Get_UpdateBoundaryConditions_Velocities(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);
        Get_UpdateHalos_Velocity(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);

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

        // Boussinesq Buoyancy Term
        //Get_BoussinesqV(MESH, V.New_Velocity, T.Pres);
		
    // 10. Calculo K3;

        Get_RK_VelocityContributions(MESH, U.K3, V.K3, W.K3);

	// 11. Calculo Nueva Velocidad

        // Intermediate Velocity Predictor U
        for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
            for(j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1]; j++){
		        for(k = 0; k < NZ; k++){
				    U.New_Velocity[LU(i,j,k,0)] = U.Pres[LU(i,j,k,0)] + RK.a43 * U.K3[LU(i,j,k,0)];
			    }
		    }
	    }

        // Intermediate Velocity Predictor V
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
			for(j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1] - 1; j++){
		    	for(k = 0; k < NZ; k++){ 
				    V.New_Velocity[LV(i,j,k,0)] = V.Pres[LV(i,j,k,0)] + RK.a43 * V.K3[LV(i,j,k,0)];
			    }
		    }
	    }

        // Intermediate Velocity Predictor W
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMW[i + Halo - Ix[Rango]][1]; j++){
		        for(k = 1; k < NZ; k++){
				    W.New_Velocity[LW(i,j,k,0)] = W.Pres[LW(i,j,k,0)] + RK.a43 * W.K3[LW(i,j,k,0)];
			    }
		    }
	    }

	// 12. Intermediate Velocities Correction
	
		Get_CorrectedVelocities(MESH, P1, U.New_Velocity, V.New_Velocity, W.New_Velocity);

	// 13. Calculo Diffusion + Convection + Boussinesq

        // Boundary Conditions Update
        Get_UpdateBoundaryConditions_Velocities(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);
        Get_UpdateHalos_Velocity(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);

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

        // Boussinesq Buoyancy Term
        //Get_BoussinesqV(MESH, V.New_Velocity, T.Pres);
		
    // 14. Calculo K4;

        Get_RK_VelocityContributions(MESH, U.K4, V.K4, W.K4);
	
    // 15. Calculo predictora final

        // Velocity Predictor U
        for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
            for(j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1]; j++){
		        for(k = 0; k < NZ; k++){
				    U.Predictor[LU(i,j,k,0)] = U.Pres[LU(i,j,k,0)] + RK.b1 * U.K1[LU(i,j,k,0)] + RK.b2 * U.K2[LU(i,j,k,0)] + RK.b3 * U.K3[LU(i,j,k,0)] + RK.b4 * U.K4[LU(i,j,k,0)];
			    }
		    }
	    }

        // Velocity Predictor V
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
			for(j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){
		    	for(k = 0; k < NZ; k++){ 
				    V.Predictor[LV(i,j,k,0)] = V.Pres[LV(i,j,k,0)] + RK.b1 * V.K1[LV(i,j,k,0)] + RK.b2 * V.K2[LV(i,j,k,0)] + RK.b3 * V.K3[LV(i,j,k,0)] + RK.b4 * V.K4[LV(i,j,k,0)];
			    }
		    }
	    }

        // Velocity Predictor W
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMW[i + Halo - Ix[Rango]][1]; j++){
		        for(k = 1; k < NZ; k++){
				    W.Predictor[LW(i,j,k,0)] = W.Pres[LW(i,j,k,0)] + RK.b1 * W.K1[LW(i,j,k,0)] + RK.b2 * W.K2[LW(i,j,k,0)] + RK.b3 * W.K3[LW(i,j,k,0)] + RK.b4 * W.K4[LW(i,j,k,0)];
			    }
		    }
	    }

		Get_UpdateBoundaryConditions_PredictorVelocities(MESH);
		
}