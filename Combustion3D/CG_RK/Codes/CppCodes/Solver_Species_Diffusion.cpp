//------------------------------------------------------------------------------------------------//
//                       CPP FILE FOR SPECIES DIFFUSION CALCULATIONS                              //
//------------------------------------------------------------------------------------------------//

// Function to calculate the binary diffusion coefficient based on Schmidt number
double Solver::Get_Diffusion_Schmidt(int n){

    return mu / (Rho * Species[n].Schmidt);

}

// Function to calculate the binary diffusion coefficient based on Lewis number
double Solver::Get_Diffusion_Lewis(int n, double T){

    return K / (Rho * JANAF_CpSpecie(T, n) * Species[n].Lewis);

}

// Function to calculate the binary diffusion coefficient based on Champan-Enskog model
double Solver::Get_Diffusion_ChampanEnskog(int sp, double Temperature, double Pressure, int i, int j, int k){
int n;
double Sum = 0.0;
double W_ab, Sigma_ab, Epsilon_ab, E_ab, Tn, OmegaD;
double D_AlphaMix;

    for (n = 0; n < N_Species; n++){
        if (n =! sp){
            W_ab = pow(1.0 / Species[sp].Wmolar + 1.0 / Species[n].Wmolar, - 1.0);
            Sigma_ab = (Species[sp].sigma + Species[n].sigma) / 2.0;
            Epsilon_ab = sqrt(Species[sp].Epsilon * Species[n].Epsilon);
            Tn = Temperature / Epsilon_ab;
            OmegaD = 1.06036 / pow(Tn, 0.15610) + 0.19300 / exp(0.47635 * Tn) + 1.03587 / exp(1.52996 * Tn) + 1.76474 / exp(3.89411 * Tn);

            D_AlphaMix = 10.1325 * (0.001858 * pow(Tn, 1.5) * pow(W_ab, -0.5)) / (Pressure * pow(Sigma_ab, 2.0) * OmegaD);

            Sum += Species[n].X[LP(i,j,k,0)] / D_AlphaMix;
        }    
    }
       
    return Sum;

}

// Function to calculate the binary diffusion coefficient based on Wilke-Lee model
double Solver::Get_Diffusion_WilkeLee(int sp, double Temperature, double Pressure, int i, int j, int k){
int n;
double Sum = 0.0;
double W_ab, Sigma_ab, Epsilon_ab, E_ab, Tn, OmegaD;
double D_AlphaMix;

    for (n = 0; n < N_Species; n++){
        if (n =! sp){
            W_ab = pow(1.0 / Species[sp].Wmolar + 1.0 / Species[n].Wmolar, - 1.0);
            Sigma_ab = (Species[sp].sigma + Species[n].sigma) / 2.0;
            Epsilon_ab = sqrt(Species[sp].Epsilon * Species[n].Epsilon);
            Tn = Temperature / Epsilon_ab;
            OmegaD = 1.06036 / pow(Tn, 0.15610) + 0.19300 / exp(0.47635 * Tn) + 1.03587 / exp(1.52996 * Tn) + 1.76474 / exp(3.89411 * Tn);

            D_AlphaMix = 10.1325 * ((0.0027 - 0.0005 * pow(W_ab, - 0.5)) * pow(Tn, 1.5) * pow(W_ab, -0.5)) / (Pressure * pow(Sigma_ab, 2.0) * OmegaD);

            Sum += Species[n].X[LP(i,j,k,0)] / D_AlphaMix;
        }    
    }
       
    return Sum;

}

// Function to calculate the diffusion coefficient based on Fick's model
void Solver::Get_Diffusion_FickModel(){
int i, j, k, n;

    for (n = 0; n < N_Species - 1; n++){
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (j = 0; j < NY; j++){
                for (k = 0; k < NZ; k++){
                    Species[n].D_AlphaMix[LP(i,j,k,0)] = (1.0 - Species[n].X[LP(i,j,k,0)]) / Get_Diffusion_ChampanEnskog(int sp, double Temperature, double Pressure, int i, int j, int k)
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
                    Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                    Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                    Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
                Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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
            Species[n].Diffusive[LP(i,j,k,0)] = Species[n].D_AlphaMix[LP(i,j,k,0)] * (1 / MESH.VolMP[LP(i,j,k,0)])*(
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