//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR CFD SOLVER CLASS UTILITIES                            //
//------------------------------------------------------------------------------------------------//

// Function to set the initial fields
void Solver::Get_InitialConditions(Mesher MESH){
int i, j, k;
double J;

	// Pressure
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				P.Pres[LP(i,j,k,0)] = 0.0;
			}
		}
	}
	
	// Velocity U
	double m = 1.7 + 0.50 * pow(Gamma_Geometry, -1.4);
    double n = 2;// + 0.3 * (1 / Gamma_Geometry - 1/3);

	for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){	
				U.Pres[LU(i,j,k,0)] = w_av * ((m+1)/m) * ((n+1)/n) * (1.0 - pow(abs(MESH.MP[LP(i,j,k,1)] - 0.50 * Ydominio) / (0.50 * Ydominio), n)) * (1.0 - pow(abs(MESH.MP[LP(i,j,k,2)] - 0.50 * Zdominio) / (0.50 * Zdominio), m));
				U.Fut[LU(i,j,k,0)] = w_av * ((m+1)/m) * ((n+1)/n) * (1.0 - pow(abs(MESH.MP[LP(i,j,k,1)] - 0.50 * Ydominio) / (0.50 * Ydominio), n)) * (1.0 - pow(abs(MESH.MP[LP(i,j,k,2)] - 0.50 * Zdominio) / (0.50 * Zdominio), m));

				U.Convective[LU(i,j,k,0)] = 0.0;
				U.Diffusive[LU(i,j,k,0)] = 0.0;	

			}
		}
	}
	
    // Velocity V
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 0; j < NY + 1; j++){
			for(k = 0; k < NZ; k++){	
				V.Pres[LV(i,j,k,0)] = 0.0;
				V.Fut[LV(i,j,k,0)] = 0.0;

				V.Convective[LV(i,j,k,0)] = 0.0;
				V.Diffusive[LV(i,j,k,0)] = 0.0;	
				V.Boussinesq[LV(i,j,k,0)] = 0.0;
				V.Boussinesq_Y[LV(i,j,k,0)] = 0.0;
			}
		}
	}

    // Velocity W
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ + 1; k++){	
				W.Pres[LW(i,j,k,0)] = 0.0;
				W.Fut[LW(i,j,k,0)] = 0.0;
				
				W.Convective[LW(i,j,k,0)] = 0.0;
				W.Diffusive[LW(i,j,k,0)] = 0.0;	
			}
		}
	}

	// Temperature T
	for(i = Ix[Rango]; i < Fx[Rango]; i++){	
		for(j = 0; j < NY; j++){
			J = j;
			for(k = 0; k < NZ; k++){	
				T.Pres[LP(i,j,k,0)] = To;
				T.Fut[LP(i,j,k,0)] = To;

				T.Convective[LP(i,j,k,0)] = 0.0;
				T.Diffusive[LP(i,j,k,0)] = 0.0;	
			}
		}
	}
	
}

// Function to calculate the mass flow entering the channel
void Solver::Get_MassFlowValue(Mesher MESH){
int i, j, k;
Qs = 0.0;

	if (Rango == 0){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){
				Qs += U.Left[LEFT(0,j,k)] * MESH.SupMU[LU(0,j,k,0)];
			}
		} 
	}

	MPI_Allreduce(&Qs, &Qs, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

}

// Function to calculate the DeltaT threshold for the diffusion
void Solver::Get_DiffusiveTimeStep(Mesher MESH){
int i, j, k;
DiffusiveDeltaT = 1000.0;

	//Comparacion 
	//a*(a>b) + b*(b>a)
	//b+=(a-b)*(a>b)

	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){	

				// CFL Diffusive U
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMU[LU(i,j,k,0)],2.0))/(mu + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMU[LU(i,j,k,0)],2.0))/(mu + 1e-10) <= DiffusiveDeltaT);
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMU[LU(i,j,k,1)],2.0))/(mu + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMU[LU(i,j,k,1)],2.0))/(mu + 1e-10) <= DiffusiveDeltaT);
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMU[LU(i,j,k,2)],2.0))/(mu + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMU[LU(i,j,k,2)],2.0))/(mu + 1e-10) <= DiffusiveDeltaT);
			
				// CFL Diffusive V
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMV[LV(i,j,k,0)],2.0))/(mu + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMV[LV(i,j,k,0)],2.0))/(mu + 1e-10) <= DiffusiveDeltaT);
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMV[LV(i,j,k,1)],2.0))/(mu + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMV[LV(i,j,k,1)],2.0))/(mu + 1e-10) <= DiffusiveDeltaT);
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMV[LV(i,j,k,2)],2.0))/(mu + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMV[LV(i,j,k,2)],2.0))/(mu + 1e-10) <= DiffusiveDeltaT);
			
				// CFL Diffusive W
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMW[LW(i,j,k,0)],2.0))/(mu + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMW[LW(i,j,k,0)],2.0))/(mu + 1e-10) <= DiffusiveDeltaT);
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMW[LW(i,j,k,1)],2.0))/(mu + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMW[LW(i,j,k,1)],2.0))/(mu + 1e-10) <= DiffusiveDeltaT);
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMW[LW(i,j,k,2)],2.0))/(mu + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMW[LW(i,j,k,2)],2.0))/(mu + 1e-10) <= DiffusiveDeltaT);
			
			}
		}
	}

	// Heat Diffusion Time Step
	for(i = Ix[Rango]; i < Fx[Rango]; i++){	
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){

				// Heat Diffusion (X Direction)
				DiffusiveDeltaT += ((0.20 * CourantFactor * Rho * Cp * pow(MESH.DeltasMP[LP(i,j,k,0)], 2.0)) / (K + 1e-10) - DiffusiveDeltaT) * ((0.20 * CourantFactor * Rho * Cp * pow(MESH.DeltasMP[LP(i,j,k,0)], 2.0)) / (K + 1e-10) <= DiffusiveDeltaT);

				// Heat Diffusion (Y Direction)
				DiffusiveDeltaT += ((0.20 * CourantFactor * Rho * Cp * pow(MESH.DeltasMP[LP(i,j,k,1)], 2.0)) / (K + 1e-10) - DiffusiveDeltaT) * ((0.20 * CourantFactor * Rho * Cp * pow(MESH.DeltasMP[LP(i,j,k,1)], 2.0)) / (K + 1e-10) <= DiffusiveDeltaT);

				// Heat Diffusion (Z Direction)
				DiffusiveDeltaT += ((0.20 * CourantFactor * Rho * Cp * pow(MESH.DeltasMP[LP(i,j,k,2)], 2.0)) / (K + 1e-10) - DiffusiveDeltaT) * ((0.20 * CourantFactor * Rho * Cp * pow(MESH.DeltasMP[LP(i,j,k,2)], 2.0)) / (K + 1e-10) <= DiffusiveDeltaT);

			}
		}
	}

	// Mass Diffusion Time Step
	for(i = Ix[Rango]; i < Fx[Rango]; i++){	
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){

				// Mass Diffusion (X Direction)
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMP[LP(i,j,k,0)],2.0))/(D_AB + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMP[LP(i,j,k,0)],2.0))/(D_AB + 1e-10) <= DiffusiveDeltaT);
				
				// Mass Diffusion (Y Direction)
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMP[LP(i,j,k,1)],2.0))/(D_AB + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMP[LP(i,j,k,1)],2.0))/(D_AB + 1e-10) <= DiffusiveDeltaT);
				
				// Mass Diffusion (Z Direction)
				DiffusiveDeltaT += ((CourantFactor * Rho*pow(MESH.DeltasMP[LP(i,j,k,2)],2.0))/(D_AB + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * Rho*pow(MESH.DeltasMP[LP(i,j,k,2)],2.0))/(D_AB + 1e-10) <= DiffusiveDeltaT);
			
			}
		}
	}
	
	MPI_Allreduce(&DiffusiveDeltaT, &DiffusiveDeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

}

//Cálculo del tiempo entre Steps
void Solver::Get_StepTime(Mesher MESH){
int i, j, k;
DeltaT = DiffusiveDeltaT;
double Tpar = 0.40;
	
	//Comparacion 
	//a*(a>b) + b*(b>a)
	//b+=(a-b)*(a>b)

	for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){	
				
				// CFL Convective U
                DeltaT += (CourantFactor * abs(MESH.DeltasMU[LU(i,j,k,0)] / (U.Pres[LU(i,j,k,0)] + 1e-10)) - DeltaT) * (CourantFactor * abs(MESH.DeltasMU[LU(i,j,k,0)] / (U.Pres[LU(i,j,k,0)] + 1e-10)) <= DeltaT);

				// CFL Convective V
				DeltaT += (CourantFactor * abs(MESH.DeltasMV[LV(i,j,k,1)] / (V.Pres[LV(i,j,k,0)] + 1e-10)) - DeltaT) * (CourantFactor * abs(MESH.DeltasMV[LV(i,j,k,1)] / (V.Pres[LV(i,j,k,0)] + 1e-10)) <= DeltaT);

				// CFL Convective W
				DeltaT += (CourantFactor * abs(MESH.DeltasMW[LW(i,j,k,2)] / (W.Pres[LW(i,j,k,0)] + 1e-10)) - DeltaT) * (CourantFactor * abs(MESH.DeltasMW[LW(i,j,k,2)] / (W.Pres[LW(i,j,k,0)] + 1e-10)) <= DeltaT);

			}	
		}
	}
	
	MPI_Allreduce(&DeltaT, &DeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

}

// Function to calculate the convective scheme
inline double Solver::ConvectiveScheme(double CoordObjetivo, double Velocity, double Coord1, double Phi1, double Coord2, double Phi2, double Coord3, double Phi3, double Coord4, double Phi4, string Esquema){

double PhiObjetivo;

double CoordD;
double PhiD;

double CoordC;
double PhiC;

double CoordU;
double PhiU;

	if (Velocity <= 0.0 || (Phi1 == 0.0 && Coord1 == 0.0)){

		CoordD = Coord2;
		PhiD = Phi2;
		CoordC = Coord3;
		PhiC = Phi3;
		CoordU = Coord4;
		PhiU = Phi4;

	}
	else if(Velocity > 0.0 || (Phi4 == 0.0 && Coord4 == 0.0)){

		CoordD = Coord3;
		PhiD = Phi3;
		CoordC = Coord2;
		PhiC = Phi2;
		CoordU = Coord1;
		PhiU = Phi1;

	}

	//Evaluacion
	double PhiF;

	if (PhiD == PhiU){
		PhiObjetivo = PhiD;
	}
	else{

		//Adimensionalizacion
		double PhiAdimC = (PhiC - PhiU)/(PhiD - PhiU + 1e-10);
		double AdimCoordC = (CoordC - CoordU)/(CoordD - CoordU + 1e-10);
		double AdimCoordE = (CoordObjetivo - CoordU)/(CoordD - CoordU + 1e-10);

        // UDS
		//PhiF = PhiAdimC;

		// SUDS
		//PhiF = (AdimCoordE / AdimCoordC) * PhiAdimC;

        // CDS
        PhiF = ((AdimCoordE - AdimCoordC) / (1 - AdimCoordC)) + ((AdimCoordE - 1) / (AdimCoordC - 1)) * PhiAdimC;

        // QUICK
		//PhiF = AdimCoordE + (((AdimCoordE*(AdimCoordE - 1.0))/(AdimCoordC*(AdimCoordC - 1.0))))*(PhiAdimC - AdimCoordC);

		//Dimensionalizacion
		PhiObjetivo = PhiU + (PhiD - PhiU)*PhiF;

	}

	return PhiObjetivo;

}

// Function to calculate the divergence of the predictors velocities (bp)
void Solver::Get_PredictorsDivergence(Mesher MESH){
int i, j, k;

    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				RHS[LA(i,j,k,0)] =  + (Rho/DeltaT)*(
								    + MESH.SupMP[LP(i,j,k,0)] * (U.Predictor[LU(i+1,j,k,0)] - U.Predictor[LU(i,j,k,0)])
								    + MESH.SupMP[LP(i,j,k,1)] * (V.Predictor[LV(i,j+1,k,0)] - V.Predictor[LV(i,j,k,0)])
								    + MESH.SupMP[LP(i,j,k,2)] * (W.Predictor[LW(i,j,k+1,0)] - W.Predictor[LW(i,j,k,0)])
								    );
			}
		}
	}

}

// Function to Check the difference between time steps
void Solver::Get_Stop(){
int i, j, k;
MaxDiffGlobal = 0.0;

	// Velocity U
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				MaxDiffGlobal += (abs((U.Fut[LU(i,j,k,0)] - U.Pres[LU(i,j,k,0)])/(U.Pres[LU(i,j,k,0)] + 1e-15)) - MaxDiffGlobal)*(abs((U.Fut[LU(i,j,k,0)] - U.Pres[LU(i,j,k,0)])/(U.Pres[LU(i,j,k,0)] + 1e-15)) >= MaxDiffGlobal);
			}
		}
	}

	// Velocity V
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 1; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				MaxDiffGlobal += (abs((V.Fut[LV(i,j,k,0)] - V.Pres[LV(i,j,k,0)])/(V.Pres[LV(i,j,k,0)] + 1e-10)) - MaxDiffGlobal)*(abs((V.Fut[LV(i,j,k,0)] - V.Pres[LV(i,j,k,0)])/(V.Pres[LV(i,j,k,0)] + 1e-10)) >= MaxDiffGlobal);
			}
		}
	}

	// Velocity W
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 1; k < NZ; k++){
				MaxDiffGlobal += (abs((W.Fut[LW(i,j,k,0)] - W.Pres[LW(i,j,k,0)])/(W.Pres[LW(i,j,k,0)] + 1e-10)) - MaxDiffGlobal)*(abs((W.Fut[LW(i,j,k,0)] - W.Pres[LW(i,j,k,0)])/(W.Pres[LW(i,j,k,0)] + 1e-10)) >= MaxDiffGlobal);
			}
		}
	}
	
	// Temperature T
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				MaxDiffGlobal += (abs((T.Fut[LP(i,j,k,0)] - T.Pres[LP(i,j,k,0)])/(T.Pres[LP(i,j,k,0)] + 1e-10)) - MaxDiffGlobal)*(abs((T.Fut[LP(i,j,k,0)] - T.Pres[LP(i,j,k,0)])/(T.Pres[LP(i,j,k,0)] + 1e-10)) >= MaxDiffGlobal);
			}
		}
	}

	// Species 0 (Water Vapour)
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				MaxDiffGlobal += (abs((Species[0].Y_Fut[LP(i,j,k,0)] - Species[0].Y_Pres[LP(i,j,k,0)])/(Species[0].Y_Pres[LP(i,j,k,0)] + 1e-10)) - MaxDiffGlobal)*(abs((Species[0].Y_Fut[LP(i,j,k,0)] - Species[0].Y_Pres[LP(i,j,k,0)])/(Species[0].Y_Pres[LP(i,j,k,0)] + 1e-10)) >= MaxDiffGlobal);
			}
		}
	}

	MPI_Allreduce(&MaxDiffGlobal, &MaxDiffGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

}

// Function to update the matrix fields
void Solver::Get_Update(){
int i, j, k;

	// Velocity U
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				U.Pres[LU(i,j,k,0)] = U.Fut[LU(i,j,k,0)];
            }
		}
	}

	// Velocity V
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 1; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				V.Pres[LV(i,j,k,0)] = V.Fut[LV(i,j,k,0)];
            }
		}
	}

	// Velocity W
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 1; k < NZ; k++){
				W.Pres[LW(i,j,k,0)] = W.Fut[LW(i,j,k,0)];
            }
		}
	}
	
	// Temperature T
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				T.Pres[LP(i,j,k,0)] = T.Fut[LP(i,j,k,0)];
				T.ContributionPast[LP(i,j,k,0)] = T.ContributionPres[LP(i,j,k,0)];
            }
		}
	}
    
	// Species 0 (Water Vapour)
	for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){
				Species[0].Y_Pres[LP(i,j,k,0)] = Species[0].Y_Fut[LP(i,j,k,0)];
				Species[0].ContributionPast[LP(i,j,k,0)] = Species[0].ContributionPres[LP(i,j,k,0)];
            }
		}
	}

}

// Function to correct the intermediate velocities
void Solver::Get_CorrectedVelocities(Mesher MESH, Parallel P1, double *U_MatrixNew, double *V_MatrixNew, double *W_MatrixNew){
int i, j, k;

	// Calculated Pressure Halo Communication
	P1.CommunicateDataLP(P.Pres, P.Pres);

	// Velocity U
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){		
				U_MatrixNew[LU(i,j,k,0)] = U_MatrixNew[LU(i,j,k,0)] - (DeltaT / Rho) * ((P.Pres[LP(i,j,k,0)] - P.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]); 
			}
		}
	}

    if (Rango == 0){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){		
				U_MatrixNew[LU(0,j,k,0)] = U.Left[LEFT(0,j,k)];
			}
		}
    }
    else if (Rango == Procesos - 1){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){		
				U_MatrixNew[LU(NX,j,k,0)] = U.Right[RIGHT(NX,j,k)];
			}
		}
    }

    // Velocity V
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 1; j < NY; j++){ 
        	for(k = 0; k < NZ; k++){		
				V_MatrixNew[LV(i,j,k,0)] = V_MatrixNew[LV(i,j,k,0)] - (DeltaT / Rho) * ((P.Pres[LP(i,j,k,0)] - P.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]);
			}
		}
	}

    // Velocity W
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 1; k < NZ; k++){		
				W_MatrixNew[LW(i,j,k,0)] = W_MatrixNew[LW(i,j,k,0)] - (DeltaT / Rho) * ((P.Pres[LP(i,j,k,0)] - P.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]);
			}
		}
	}

}

// Function to calculate the new corrected velocities
void Solver::Get_Velocities(Mesher MESH, Parallel P1){
int i, j, k;

	// Calculated Pressure Halo Communication
	P1.CommunicateDataLP(P.Pres, P.Pres);

	// Velocity U
    for(i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){		
				U.Fut[LU(i,j,k,0)] = U.Predictor[LU(i,j,k,0)] - (DeltaT / Rho) * ((P.Pres[LP(i,j,k,0)] - P.Pres[LP(i-1,j,k,0)]) / MESH.DeltasMU[LU(i,j,k,0)]); 
			}
		}
	}

    if (Rango == 0){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){		
				U.Fut[LU(0,j,k,0)] = U.Left[LEFT(0,j,k)];
			}
		}
    }
    
	if (Rango == Procesos - 1){
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){		
				U.Fut[LU(NX,j,k,0)] = U.Right[RIGHT(NX,j,k)];
			}
		}
    }

    // Velocity V
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 1; j < NY; j++){
        	for(k = 0; k < NZ; k++){ 		
				V.Fut[LV(i,j,k,0)] = V.Predictor[LV(i,j,k,0)] - (DeltaT / Rho) * ((P.Pres[LP(i,j,k,0)] - P.Pres[LP(i,j-1,k,0)]) / MESH.DeltasMV[LV(i,j,k,1)]);
			}
		}
	}

    // Velocity W
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
		    for(k = 1; k < NZ; k++){		
				W.Fut[LW(i,j,k,0)] = W.Predictor[LW(i,j,k,0)] - (DeltaT / Rho) * ((P.Pres[LP(i,j,k,0)] - P.Pres[LP(i,j,k-1,0)]) / MESH.DeltasMW[LW(i,j,k,2)]);
			}
		}
	}

	Get_UpdateBoundaryConditions_Velocities(U.Fut, V.Fut, W.Fut);
	Get_UpdateHalos_Velocity(U.Fut, V.Fut, W.Fut);

}