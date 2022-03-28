//------------------------------------------------------------------------------------------------//
//                       CPP FILE FOR SPECIES DIFFUSION CALCULATIONS                              //
//------------------------------------------------------------------------------------------------//

// Function to calculate species binary diffusion coefficients
void Solver::Get_Species_DiffusionCoefficients(){
int i, j, k, n;

    for (n = 0; n < N_Species - 1; n++){
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    Species[n].D_ab[LP(i,j,k,0)] = D_AB;
                }
            }
        }
    }
    
    
}

// Function to calculate the species diffusion
void Solver::Get_Species_Diffusion(Mesher MESH){
int i, j, k, n;

    for (n = 0; n < N_Species - 1; n++){
      
        // Center
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = 1; j < NY - 1; j++){
                for (k = 1; k < NZ - 1; k++){
                    Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
                }
            }
        }

        // Bottom
        j = 0;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (k = 1; k < NZ - 1; k++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }
        }

        // Top
        j = NY - 1;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (k = 1; k < NZ - 1; k++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Top[TOP(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }
        }

        // Here
        k = 0;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = 1; j < NY - 1; j++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }
        }

        // There
        k = NZ - 1;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = 1; j < NY - 1; j++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].There[THERE(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }
        }

        // Bottom Here Corner
        j = 0;
        k = 0;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }

        // Bottom There Corner
        j = 0;
        k = NZ - 1;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].There[THERE(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }

        // Top Here Corner
        j = NY - 1;
        k = 0;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Top[TOP(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }


        // Top There Corner
        i = NY - 1;
        k = NZ - 1;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Top[TOP(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].There[THERE(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
        }

        if (Rango == 0){
            i = 0;

            // Center
            for (j = 1; j < NY - 1; j++){
                for (k = 1; k < NZ - 1; k++){
                    Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
                }
            }
    
            // Bottom
            j = 0;
            for (k = 1; k < NZ - 1; k++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }

            // Top
            j = NY - 1;
            for (k = 1; k < NZ - 1; k++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Top[TOP(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }

            // Here
            k = 0;
            for (j = 1; j < NY - 1; j++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }

            // There
            k = NZ - 1;
            for (j = 1; j < NY - 1; j++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].There[THERE(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }

            // Bottom Here Corner
            j = 0;
            k = 0;
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
    

            // Bottom There Corner
            j = 0;
            k = NZ - 1;
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].There[THERE(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
    

            // Top Here Corner
            j = NY - 1;
            k = 0;
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Top[TOP(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
    


            // Top There Corner
            i = NY - 1;
            k = NZ - 1;
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i+1,j,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Left[LEFT(i,j,k)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Top[TOP(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].There[THERE(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );

        }
        else if (Rango == Procesos - 1){
            i = NX - 1;

            // Center
            for (j = 1; j < NY - 1; j++){
                for (k = 1; k < NZ - 1; k++){
                    Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Right[RIGHT(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
                }
            }

            // Bottom
            j = 0;
            for (k = 1; k < NZ - 1; k++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Right[RIGHT(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }

            // Top
            j = NY - 1;
            for (k = 1; k < NZ - 1; k++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Right[RIGHT(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Top[TOP(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }

            // Here
            k = 0;
            for (j = 1; j < NY - 1; j++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Right[RIGHT(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }

            // There
            k = NZ - 1;
            for (j = 1; j < NY - 1; j++){
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Right[RIGHT(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].There[THERE(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
            }

            // Bottom Here Corner
            j = 0;
            k = 0;
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Right[RIGHT(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
    

            // Bottom There Corner
            j = 0;
            k = NZ - 1;
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Right[RIGHT(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j+1,k,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Bottom[BOTTOM(i,j,k)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].There[THERE(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
    

            // Top Here Corner
            j = NY - 1;
            k = 0;
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Right[RIGHT(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Top[TOP(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k+1,0)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Here[HERE(i,j,k)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
    


            // Top There Corner
            i = NY - 1;
            k = NZ - 1;
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_ab[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
                                                   + MESH.SupMP[LP(i,j,k,0)] * (Species[n].Right[RIGHT(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMU[LU(i+1,j,k,0)]
                                                   - MESH.SupMP[LP(i,j,k,0)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]
                                                   + MESH.SupMP[LP(i,j,k,1)] * (Species[n].Top[TOP(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMV[LV(i,j+1,k,1)]
                                                   - MESH.SupMP[LP(i,j,k,1)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]
                                                   + MESH.SupMP[LP(i,j,k,2)] * (Species[n].There[THERE(i,j,k)] - Species[n].Y_Pres[LP(i,j,k,0)]) / MESH.DeltasMW[LW(i,j,k+1,2)]
                                                   - MESH.SupMP[LP(i,j,k,2)] * (Species[n].Y_Pres[LP(i,j,k,0)] - Species[n].Y_Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]
										           );
    

        }
    }
    
}