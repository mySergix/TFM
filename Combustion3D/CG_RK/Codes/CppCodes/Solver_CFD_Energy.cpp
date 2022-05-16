//------------------------------------------------------------------------------------------------//
//                         CPP FILE FOR N-S ENERGY EQUATION CALCULATIONS                          //
//------------------------------------------------------------------------------------------------//

// Function to calculaHetTe diffusive term of energy equation
void Solver::Get_DiffusionEnergy(Mesher MESH){
int i, j, k;

    // Center
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){
            for (k = 1; k < NZ - 1; k++){

                // Core
                T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

                // Internal Left
				if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                    T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (Twalls_IntLeft - T.Pres[LP(i,j,k,0)]) / (0.50 * MESH.DeltasMU[LU(i+1,j,k,0)])
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / (0.50 * MESH.DeltasMU[LU(i,j,k,0)])
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
                }
                // Internal Right
				else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                    T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / (0.50 * MESH.DeltasMU[LU(i+1,j,k,0)])
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - Twalls_IntRight) / (0.50 * MESH.DeltasMU[LU(i,j,k,0)])
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
                }

            }
        }
    }

    // Bottom
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
        for (k = 1; k < NZ - 1; k++){

            // Core
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
                                        
            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (Twalls_IntLeft - T.Pres[LP(i,j,k,0)]) / (0.50 * MESH.DeltasMU[LU(i+1,j,k,0)])
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / (0.50 * MESH.DeltasMU[LU(i+1,j,k,0)])
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - Twalls_IntRight) / (0.50 * MESH.DeltasMU[LU(i,j,k,0)])
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }

        }
    }

    // Top
    j = NY - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 1; k < NZ - 1; k++){
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Top[TOP(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }
    }

    // Here
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){

            // Core
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (Twalls_IntLeft - T.Pres[LP(i,j,k,0)]) / (0.50 * MESH.DeltasMU[LU(i+1,j,k,0)])
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - Twalls_IntRight) / (0.50 * MESH.DeltasMU[LU(i,j,k,0)])
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }

        }
    }

    // There
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){

            // Core
            T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (Twalls_IntLeft - T.Pres[LP(i,j,k,0)]) / (0.50 * MESH.DeltasMU[LU(i+1,j,k,0)])
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - Twalls_IntRight) / (0.50 * MESH.DeltasMU[LU(i,j,k,0)])
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }

        }
    }

    // Bottom Here Corner 
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];

            // Core
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (Twalls_IntLeft - T.Pres[LP(i,j,k,0)]) / (0.50 * MESH.DeltasMU[LU(i+1,j,k,0)])
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - Twalls_IntRight) / (0.50 * MESH.DeltasMU[LU(i,j,k,0)])
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }

    }

    // Bottom There Corner
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];

            // Core
            T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (Twalls_IntLeft - T.Pres[LP(i,j,k,0)]) / (0.50 * MESH.DeltasMU[LU(i+1,j,k,0)])
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - Twalls_IntRight) / (0.50 * MESH.DeltasMU[LU(i,j,k,0)])
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }

    }

    // Top Here Corner
    j = NY - 1;
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){

            // Core
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Top[TOP(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

    }

    // Top There Corner
    j = NY - 1;
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){

            // Core
            T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Top[TOP(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );

    }

    if (Rango == 0){
        i = 0;

        // Center
        for (j = 1; j < NY - 1; j++){
            for (k = 1; k < NZ - 1; k++){
                T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }
        }

        // Bottom
        j = 0;
        for (k = 1; k < NZ - 1; k++){
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }

        // Top
        j = NY - 1;
        for (k = 1; k < NZ - 1; k++){
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Top[TOP(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }

        // Here
        k = 0;
        for (j = 1; j < NY - 1; j++){
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }

        // There
        k = NZ - 1;
        for (j = 1; j < NY - 1; j++){
            T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }

        // Bottom Here Corner
        j = 0;
        k = 0;
        
                T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        

        // Bottom There Corner
        j = 0;
        k = NZ - 1;
        
            T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        


        // Top Here Corner
        j = NY - 1;
        k = 0;
        
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Top[TOP(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        

        // Top There Corner
        j = NY - 1;
        k = NZ - 1;
        
            T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + 0.50 * (K_Thermal[LP(i+1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Top[TOP(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        

    }
    else if (Rango == Procesos - 1){
        i = NX - 1;

        // Center
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){
            for (k = 1; k < NZ - 1; k++){
                T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Right[RIGHT(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
            }
        }

        // Bottom
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
        for (k = 1; k < NZ - 1; k++){
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Right[RIGHT(I,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }

        // Top
        j = NY - 1;
        for (k = 1; k < NZ - 1; k++){
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i+1,j,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Top[TOP(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }

        // Here
        k = 0;
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Right[RIGHT(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }

        // There
        k = NZ - 1;
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){
            T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Right[RIGHT(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        }

        // Bottom Here Corner
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
        k = 0;
        
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Right[RIGHT(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        

        // Bottom There Corner
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
        k = NZ - 1;
        
            T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Right[RIGHT(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + 0.50 * (K_Thermal[LP(i,j+1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j+1,k,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        


        // Top Here Corner
        j = NY - 1;
        k = 0;
        
            T.Diffusive[LP(i,j,k,0)] = (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Right[RIGHT(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Top[TOP(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + 0.50 * (K_Thermal[LP(i,j,k+1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k+1,0)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        

        // Top There Corner
        j = NY - 1;
        k = NZ - 1;
            T.Diffusive[LP(i,j,k,0)] =  (1.0 / (Rho * MESH.VolMP[LP(i,j,k,0)]))*(
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,0)] * (T.Right[RIGHT(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                         - 0.50 * (K_Thermal[LP(i-1,j,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,0)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,1)] * (T.Top[TOP(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                         - 0.50 * (K_Thermal[LP(i,j-1,k,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,1)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                         + K_Thermal[LP(i,j,k,0)] * MESH.SupMP[LP(i,j,k,2)] * (T.There[THERE(i,j,k)] - T.Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                         - 0.50 * (K_Thermal[LP(i,j,k-1,0)] + K_Thermal[LP(i,j,k,0)]) * MESH.SupMP[LP(i,j,k,2)] * (T.Pres[LP(i,j,k,0)] - T.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										 );
        

    }

}

// Function to calculate tTe convective term of energy equation
void Solver::Get_ConvectionEnergy(Mesher MESH){
int i, j, k;
double Tw, Te, Ts, Tn, Th, Tt;

    // Center
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){
            for (k = 1; k < NZ - 1; k++){

                // Core
                Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                // Internal Left
				if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                    Te = Twalls_IntLeft;
                }
                // Internal Right
				else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                    Tw = Twalls_IntRight;
                }

                T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

            }
        }
    }

    // Bottom
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
        for (k = 1; k < NZ - 1; k++){

            // Core
            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = T.Bottom[BOTTOM(i,j,k)];
            Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                Te = Twalls_IntLeft;
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                Tw = Twalls_IntRight;
            }

            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Bottom[BOTTOM(i,j,k)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

        }
    }

    // Top
    j = NY - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 1; k < NZ - 1; k++){

            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Tn = T.Top[TOP(i,j,k)];

            Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Top[TOP(i,j,k)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

        }
    }

    // Here
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){
            
            // Core
            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Th = T.Here[HERE(i,j,k)];
            Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                Te = Twalls_IntLeft;
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                Tw = Twalls_IntRight;
            }

            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Here[HERE(i,j,k)])
											 );

        }
    }

    // There
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){

            // Core
            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Tt = T.There[THERE(i,j,k)];

            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                Te = Twalls_IntLeft;
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                Tw = Twalls_IntRight;
            }

            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.There[THERE(i,j,k)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

        }
    }

    // Bottom Here Corner
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];

            // Core
            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = T.Bottom[BOTTOM(i,j,k)];
            Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Th = T.Here[HERE(i,j,k)];
            Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                Te = Twalls_IntLeft;
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                Tw = Twalls_IntRight;
            }

            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Here[HERE(i,j,k)])
											 );

    }

    // Bottom There Corner  
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];

            // Core
            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = T.Bottom[BOTTOM(i,j,k)];
            Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Tt = T.There[THERE(i,j,k)];

            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                Te = Twalls_IntLeft;
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                Tw = Twalls_IntRight;
            }

            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.There[THERE(i,j,k)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

    }

    // Top Here Corner
    j = NY - 1;
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){

            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Tn = T.Top[TOP(i,j,k)];
            
            Th = T.Here[HERE(i,j,k)];
            Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Here[HERE(i,j,k)])
											 );

    }

    // Top There Corner
    j = NY - 1;
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){

            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Tn = T.Top[TOP(i,j,k)];
            
            Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Tt = T.There[THERE(i,j,k)];

            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.There[THERE(i,j,k)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

    }

    if (Rango == 0){
        i = 0;

        // Center
        for (j = 1; j < NY - 1; j++){
            for (k = 1; k < NZ - 1; k++){

                Tw = T.Left[LEFT(i,j,k)];
                Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Left[LEFT(i,j,k)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

            }
        }

        // Bottom
        j = 0;
        for (k = 1; k < NZ - 1; k++){

            Tw = T.Left[LEFT(i,j,k)];
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = T.Bottom[BOTTOM(i,j,k)];
            Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Left[LEFT(i,j,k)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Bottom[BOTTOM(i,j,k)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

        }

        // Top
        j = NY - 1;
        for (k = 1; k < NZ - 1; k++){

            Tw = T.Left[LEFT(i,j,k)];
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Tn = T.Top[TOP(i,j,k)];

            Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Left[LEFT(i,j,k)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Top[TOP(i,j,k)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

        }

        // Here
        k = 0;
        for (j = 1; j < NY - 1; j++){

            Tw = T.Left[LEFT(i,j,k)];
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Th = T.Here[HERE(i,j,k)];
            Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Left[LEFT(i,j,k)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Here[HERE(i,j,k)])
											 );

        }

        // There
        k = NZ - 1;
        for (j = 1; j < NY - 1; j++){

            Tw = T.Left[LEFT(i,j,k)];
            Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Tt = T.There[THERE(i,j,k)];

            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Left[LEFT(i,j,k)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.There[THERE(i,j,k)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

        }

        // Bottom Here Corner
        j = 0;
        k = 0;

        Tw = T.Left[LEFT(i,j,k)];
        Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
        Ts = T.Bottom[BOTTOM(i,j,k)];
        Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
        Th = T.Here[HERE(i,j,k)];
        Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
        T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Here[HERE(i,j,k)])
											 );


        // Bottom There Corner
        j = 0;
        k = NZ - 1;
    
        Tw = T.Left[LEFT(i,j,k)];
        Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
        Ts = T.Bottom[BOTTOM(i,j,k)];
        Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
        Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
        Tt = T.There[THERE(i,j,k)];

        T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.There[THERE(i,j,k)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

    

        // Top Here Corner
        j = NY - 1;
        k = 0;
    
        Tw = T.Left[LEFT(i,j,k)];
        Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
        Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
        Tn = T.Top[TOP(i,j,k)];
            
        Th = T.Here[HERE(i,j,k)];
        Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
        T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Here[HERE(i,j,k)])
											 );

    

        // Top There Corner
        k = NZ - 1;
    
        Tw = T.Left[LEFT(i,j,k)];
        Te = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], T.Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
        Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
        Tn = T.Top[TOP(i,j,k)];
            
        Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
        Tt = T.There[THERE(i,j,k)];

        T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.There[THERE(i,j,k)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

    

    }
    else if (Rango == Procesos - 1){
        i = NX - 1;

        // Center
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){
            for (k = 1; k < NZ - 1; k++){

                Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Te = T.Right[RIGHT(i,j,k)];

                Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Right[RIGHT(i,j,k)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

            }
        }

        // Bottom
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
        for (k = 1; k < NZ - 1; k++){

            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = T.Right[RIGHT(i,j,k)];

            Ts = T.Bottom[BOTTOM(i,j,k)];
            Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Right[RIGHT(i,j,k)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Bottom[BOTTOM(i,j,k)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

        }

        // Top
        j = NY - 1;
        for (k = 1; k < NZ - 1; k++){

            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = T.Right[RIGHT(i,j,k)];

            Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Tn = T.Top[TOP(i,j,k)];

            Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Right[RIGHT(i,j,k)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Top[TOP(i,j,k)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

        }

        // Here
        k = 0;
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){

            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = T.Right[RIGHT(i,j,k)];

            Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Th = T.Here[HERE(i,j,k)];
            Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Right[RIGHT(i,j,k)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Here[HERE(i,j,k)])
											 );

        }

        // There
        k = NZ - 1;
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){

            Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Te = T.Right[RIGHT(i,j,k)];

            Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Tt = T.There[THERE(i,j,k)];

            T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Right[RIGHT(i,j,k)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.There[THERE(i,j,k)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

        }

        // Bottom Here Corner
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
        k = 0;
    
        Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
        Te = T.Right[RIGHT(i,j,k)];

        Ts = T.Bottom[BOTTOM(i,j,k)];
        Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
        Th = T.Here[HERE(i,j,k)];
        Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
        T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Here[HERE(i,j,k)])
											 );

    

        // Bottom There Corner
        j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
        k = NZ - 1;
    
        Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
        Te = T.Right[RIGHT(i,j,k)];

        Ts = T.Bottom[BOTTOM(i,j,k)];
        Tn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], T.Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
        Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
        Tt = T.There[THERE(i,j,k)];

        T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.There[THERE(i,j,k)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

        // Top Here Corner
        j = NY - 1;
        k = 0;
    
        Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
        Te = T.Right[RIGHT(i,j,k)];

        Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
        Tn = T.Top[TOP(i,j,k)];
            
        Th = T.Here[HERE(i,j,k)];
        Tt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], T.Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
        T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Tt - Th * W.Here[HERE(i,j,k)])
											 );

        // Top There Corner
        j = NY - 1;
        k = NZ - 1;
    
        Tw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], T.Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], T.Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], T.Pres[LP(i+1,j,k,0)], EsquemaLargo);
        Te = T.Right[RIGHT(i,j,k)];

        Ts = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], T.Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], T.Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], T.Pres[LP(i,j+1,k,0)], EsquemaLargo);
        Tn = T.Top[TOP(i,j,k)];
            
        Th = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], T.Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], T.Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], T.Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], T.Pres[LP(i,j,k+1,0)], EsquemaLargo);
        Tt = T.There[THERE(i,j,k)];

        T.Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Te - Tw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Tn - Ts * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.There[THERE(i,j,k)] * Tt - Th * W.Pres[LW(i,j,k,0)])
											 );

    }

}

// Function to calculaHetTe enthalpy variation due to chemical reactions
void Solver::Get_ReactionsEnergy(Mesher MESH){
int i, j, k, n;
 
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
            for (k = 0; k < NZ; k++){
                T.Reactive[LP(i,j,k,0)] = 0.0;
                for (n = 0; n < N_Species - 1; n++){
                    T.Reactive[LP(i,j,k,0)] += Species[n].wk[LP(i,j,k,0)] * JANAF_AbsEnthalpy_Specie(n, T.Pres[LP(i,j,k,0)]);
                }
            }
        }
    }

}

// Function to calculate tTe contributions of tTe energy equation
void Solver::Get_EnergyContributions(Mesher MESH){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
            for (k = 0; k < NZ; k++){
                T.ContributionPres[LP(i,j,k,0)] = (T.Diffusive[LP(i,j,k,0)] + T.Reactive[LP(i,j,k,0)]) / 1000.0 - T.Convective[LP(i,j,k,0)];
            }
        }
    }

}

// Function to calculate tTe new temperatures field
void Solver::Get_NewTemperatures(Mesher MESH){
int i, j, k;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
            for (k = 0; k < NZ; k++){
                T.Fut[LP(i,j,k,0)] = T.Pres[LP(i,j,k,0)] + DeltaT * (1.50 * T.ContributionPres[LP(i,j,k,0)] - 0.50 * T.ContributionPast[LP(i,j,k,0)]);
            }
        }
    }

}