//------------------------------------------------------------------------------------------------//
//                CPP FILE FOR N-S MOMENTUM EQUATION DIFFUSION CALCULATIONS                       //
//------------------------------------------------------------------------------------------------//

// Function to calculate the diffusion term of U Velocity
void Solver::Get_DiffusionU(Mesher MESH, double *U_Matrix){
int i, j, k;

    // Centro
	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = NY_ColumnMU[i + Halo - Ix[Rango]][0] + 1; j < NY_ColumnMU[i + Halo - Ix[Rango]][1] - 1; j++){
		    for(k = 1; k < NZ - 1; k++){	

				// Core
				U.Diffusive[LU(i,j,k,0)] =  (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

				// Internal Left
				if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
					U.Diffusive[LU(i,j,k,0)] =  (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.I_Left[ILEFT(0,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
				}
				// Internal Right
				else if (i == NX_1 + NX_2 + 1 && j < NY - (NY_3 + NY_4)){
					U.Diffusive[LU(i,j,k,0)] =  (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U.I_Right[IRIGHT(0,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
				}

			}
		}
	}
	
	
    // Top
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		j = NY_ColumnMU[i + Halo - Ix[Rango]][1] - 1;
		for(k = 1; k < NZ - 1; k++){
			U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Top[TOP(i,j,k)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }
    }
	
	// Bottom
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		j = NY_ColumnMU[i + Halo - Ix[Rango]][0];
		for(k = 1; k < NZ - 1; k++){

			// Core
			U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

			// Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
				U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.I_Left[ILEFT(0,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
			// Internal Right
			else if (i == NX_1 + NX_2 + 1 && j < NY - (NY_3 + NY_4)){
				U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U.I_Right[IRIGHT(0,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
            
        }
    }

	// Here
    k = 0;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = NY_ColumnMU[i + Halo - Ix[Rango]][0] + 1; j < NY_ColumnMU[i + Halo - Ix[Rango]][1] - 1; j++){

			// Core
			U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			
			// Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
				U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.I_Left[ILEFT(0,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
			// Internal Right
			else if (i == NX_1 + NX_2 + 1 && j < NY - (NY_3 + NY_4)){
				U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U.I_Right[IRIGHT(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}

        }
    }

	// There
    k = NZ - 1;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = NY_ColumnMU[i + Halo - Ix[Rango]][0] + 1; j < NY_ColumnMU[i + Halo - Ix[Rango]][1] - 1; j++){

			// Core
			U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.There[THERE(i,j,k)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			
			// Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
				U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.I_Left[ILEFT(0,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.There[THERE(i,j,k)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
			// Internal Right
			else if (i == NX_1 + NX_2 + 1 && j < NY - (NY_3 + NY_4)){
				U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U.I_Right[IRIGHT(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.There[THERE(i,j,k)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
            
        }
    }

    // Bottom Here Corner
    k = 0;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		j = NY_ColumnMU[i + Halo - Ix[Rango]][0];

		// Core
		U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		
		// Internal Left
		if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
			U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.I_Left[ILEFT(0,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}
		// Internal Right
		else if (i == NX_1 + NX_2 + 1 && j < NY - (NY_3 + NY_4)){
			U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U.I_Right[IRIGHT(0,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}
  
    }

    // Bottom There Corner
    k = NZ - 1;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		j = NY_ColumnMU[i + Halo - Ix[Rango]][0];

		// Core
		U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.There[THERE(i,j,k)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		
		// Internal Left
		if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
			U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U.I_Left[ILEFT(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.There[THERE(i,j,k)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}
		// Internal Right
		else if (i == NX_1 + NX_2 + 1 && j < NY - (NY_3 + NY_4)){
			U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U.I_Right[IRIGHT(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j+1,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.There[THERE(i,j,k)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}		
        
    }

    // Top Here Corner
    k = 0;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		j = NY_ColumnMU[i + Halo - Ix[Rango]][1] - 1;

        U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Top[TOP(i,j,k)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k+1,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
    }

    // Top There Corner
    k = NZ - 1;
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		j = NY_ColumnMU[i + Halo - Ix[Rango]][1] - 1;

        U.Diffusive[LU(i,j,k,0)] = (mu/(Rho*MESH.VolMU[LU(i,j,k,0)]))*(
                                         + MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i+1,j,k,0)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,0)]
                                         - MESH.SupMU[LU(i,j,k,0)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i-1,j,k,0)]) / MESH.DeltasMP[LP(i-1,j,k,0)]
										 + MESH.SupMU[LU(i,j,k,1)] * (U.Top[TOP(i,j,k)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMU[LU(i,j,k,1)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMU[LU(i,j,k,2)] * (U.There[THERE(i,j,k)] - U_Matrix[LU(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
										 - MESH.SupMU[LU(i,j,k,2)] * (U_Matrix[LU(i,j,k,0)] - U_Matrix[LU(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
    }
	

    if (Rango == 0){
		i = 0;
        for(j = NY_ColumnMU[i + Halo - Ix[Rango]][0]; j < NY_ColumnMU[i + Halo - Ix[Rango]][1]; j++){
		    for(k = 0; k < NZ; k++){	
				U.Diffusive[LU(0,j,k,0)] = 0.0;
			}
		}
    }
    else if (Rango == Procesos - 1){
		i = NX;
        for(j = NY_ColumnMU[i + Halo - Ix[Rango]][0]; j < NY_ColumnMU[i + Halo - Ix[Rango]][1]; j++){
		    for(k = 0; k < NZ; k++){	
				U.Diffusive[LU(NX,j,k,0)] = 0.0;
			}
		}
    }

}

// Function to calculate the diffusion term of V Velocity
void Solver::Get_DiffusionV(Mesher MESH, double *V_Matrix){
int i, j, k;

    // Center
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){
			for(k = 1; k < NZ - 1; k++){
			
				// Core 
				V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i+1,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k+1,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

				// Internal Left
				if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
					V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.I_Left[ILEFT(0,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k+1,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
				}
				// Internal Right
				else if (i == NX_1 + NX_2 + 1 && j < NY - (NY_3 + NY_4)){
					V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i+1,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V.I_Right[IRIGHT(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k+1,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
				}

			}
		}
	}

    // Here
    k = 0;
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){

			// Core
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i+1,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k+1,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

			// Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
				V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.I_Left[ILEFT(0,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k+1,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
			// Internal Right
			else if (i == NX_1 + NX_2 + 1 && j < NY - (NY_3 + NY_4)){
				V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i+1,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V.I_Right[IRIGHT(0,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k+1,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}

		}
	}

    // There
    k = NZ - 1;
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){

			// Core
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i+1,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.There[THERE(i,j,k)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

			// Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
				V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.I_Left[ILEFT(i+1,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.There[THERE(i,j,k)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
			// Internal Right
			else if (i == NX_1 + NX_2 + 1 && j < NY - (NY_3 + NY_4)){
				V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i+1,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V.I_Right[IRIGHT(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.There[THERE(i,j,k)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}

		}
	}

    if (Rango == 0){

        i = 0;
		for(k = 1; k < NZ - 1; k++){
			for(j = 1; j < NY; j++){
				V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i+1,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k+1,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
		}

        // Here
        k = 0;
		for(j = 1; j < NY; j++){
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i+1,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k+1,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}

        // There
        k = NZ - 1;
		for(j = 1; j < NY; j++){
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i+1,j,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.There[THERE(i,j,k)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}

    }
    else if (Rango == Procesos - 1){

        i = NX - 1;
		for(k = 1; k < NZ - 1; k++){
			for(j = NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){
				V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Right[RIGHT(i,j,k)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k+1,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
			}
		}

        // Here
        k = 0;
		for(j = NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Right[RIGHT(i,j,k)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k+1,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}

        // There
        k = NZ - 1;
		for(j = NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){
			V.Diffusive[LV(i,j,k,0)] = (mu/(Rho*MESH.VolMV[LV(i,j,k,0)]))*(
										 + MESH.SupMV[LV(i,j,k,0)] * (V.Right[RIGHT(i,j,k)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMV[LV(i,j,k,0)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]	
										 + MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j+1,k,0)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,1)]
                                         - MESH.SupMV[LV(i,j,k,1)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j-1,k,0)]) / MESH.DeltasMP[LP(i,j-1,k,1)]
										 + MESH.SupMV[LV(i,j,k,2)] * (V.There[THERE(i,j,k)] - V_Matrix[LV(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - MESH.SupMV[LV(i,j,k,2)] * (V_Matrix[LV(i,j,k,0)] - V_Matrix[LV(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
		}

	}

}

// Function to calculate the diffusion term of W Velocity
void Solver::Get_DiffusionW(Mesher MESH, double *W_Matrix){
int i, j, k;

    // Center
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = NY_ColumnMW[i + Halo - Ix[Rango]][0] + 1; j < NY_ColumnMW[i + Halo - Ix[Rango]][1] - 1; j++){
			for(k = 1; k < NZ; k++){

				// Core
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i+1,j,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j+1,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );

				// Internal Left
				if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
					W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.I_Left[ILEFT(0,j,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j+1,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
				}
				// Internal Right
				else if (i == NX_1 + NX_2 + 1 && j < NY - (NY_3 + NY_4)){
					W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i+1,j,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W.I_Right[IRIGHT(0,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j+1,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
				}

			}
		}
    }

    // Bottom
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		j = NY_ColumnMW[i + Halo - Ix[Rango]][0];
		for(k = 1; k < NZ; k++){

			// Core
			W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i+1,j,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j+1,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );

			// Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.I_Left[ILEFT(0,j,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j+1,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
			}
			// Internal Right
			else if (i == NX_1 + NX_2 + 1 && j < NY - (NY_3 + NY_4)){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i+1,j,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W.I_Right[IRIGHT(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j+1,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
			}

		}
    }

    // Top
    j = NY - 1;
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(k = 1; k < NZ; k++){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i+1,j,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Top[TOP(i,j,k)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
		}
    }

    if (Rango == 0){

        i = 0;
		for(k = 1; k < NZ; k++){
			for(j = 1; j < NY - 1; j++){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i+1,j,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j+1,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
			}
		}

        // Bottom
        j = 0;
		for(k = 1; k < NZ; k++){
			W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i+1,j,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j+1,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
		}
    
        // Top
        j = NY - 1;
		for(k = 1; k < NZ; k++){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i+1,j,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Top[TOP(i,j,k)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
		}

    }
    else if(Rango == Procesos - 1){

        i = NX - 1;
		for(k = 1; k < NZ; k++){
			for(j = NY_ColumnMW[i + Halo - Ix[Rango]][0] + 1; j < NY_ColumnMW[i + Halo - Ix[Rango]][1] - 1; j++){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Right[RIGHT(i,j,k)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j+1,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
			}
		}

        // Bottom
        j = NY_ColumnMW[i + Halo - Ix[Rango]][0];
		for(k = 1; k < NZ; k++){
			W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Right[RIGHT(i,j,k)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j+1,k,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
		}

        // Top
        j = NY_ColumnMW[i + Halo - Ix[Rango]][1] - 1;
		for(k = 1; k < NZ; k++){
				W.Diffusive[LW(i,j,k,0)] = (mu/(Rho*MESH.VolMW[LW(i,j,k,0)]))*(
										 + MESH.SupMW[LW(i,j,k,0)] * (W.Right[RIGHT(i,j,k)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - MESH.SupMW[LW(i,j,k,0)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
										 + MESH.SupMW[LW(i,j,k,1)] * (W.Top[TOP(i,j,k)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - MESH.SupMW[LW(i,j,k,1)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
										 + MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k+1,0)] - W_Matrix[LW(i,j,k,0)]) / MESH.DeltasMP[LP(i,j,k,2)]
										 - MESH.SupMW[LW(i,j,k,2)] * (W_Matrix[LW(i,j,k,0)] - W_Matrix[LW(i,j,k-1,0)]) / MESH.DeltasMP[LP(i,j,k-1,2)]
										 );
		}

    }

}
