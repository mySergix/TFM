//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR SOLVER EXECUTION                                      //
//------------------------------------------------------------------------------------------------//

// Function to run the solver of the simulation
void Solver::RunSolver(Memory M1, ReadData R1, Parallel P1, Mesher MESH, PostProcessing POST1){	
int i, j, k;

	double Time = 0.0;
	int Step = 0;

	char FileName_1[300]; 

	MaxDiffGlobal = 2.0*ConvergenciaGlobal;

	// Runge Kutta 3rd Order Coefficients
	Get_RK_Coefficients(RK3);

	// Memory Allocation for Navier Stokes Equations
	Allocate_VelocitiesMemory(M1);
    Allocate_PressureMemory(M1);
    Allocate_PoissonCoeffsMemory(M1);
	Allocate_EnergyMemory(M1);
	if (Rango == 0){ Allocate_GlobalMemory(M1); }

    // Initial Simulation settings
    Get_StaticBoundaryConditions_Velocities(MESH);
	Get_StaticHalos_Velocity();

	if (Problema == 2){
		Get_StaticBoundaryConditions_Temperatures(MESH);
		Get_StaticHalos_Temperatures();
	}

	// Initial Calculations
    Get_PoissonCoefficients(MESH);	
    
	// Communication to global matrix at step 0
	P1.SendMatrixToZeroMP(P.Pres, Global.P);
	P1.SendMatrixToZeroMP(T.Pres, Global.T);
	P1.SendMatrixToZeroMU(U.Fut, Global.U);
	P1.SendMatrixToZeroMV(V.Fut, Global.V);
	P1.SendMatrixToZeroMW(W.Fut, Global.W);

	// Print of Step 0 .VTK files
	if(Rango == 0){
		sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("DrivenCavity/", "Presion", FileName_1, MESH, Global.P);

		if (Problema == 2){
			sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
			POST1.VTK_GlobalScalar3D("DrivenCavity/", "Temperatura", FileName_1, MESH, Global.T);
		}
		
		sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
		POST1.VTK_GlobalVectorial3D("DrivenCavity/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);
	}

	while(Step < 10){
		//MaxDiffGlobal >= ConvergenciaGlobal
		// New Step
		Step++;
		
		// Step Time Calculation
		Get_StepTime(MESH, P1);
		Time += DeltaT;
        
		// RK3 Integration (Velocity + Temperature)
		Get_RK3_Integration(MESH, P1);

		// Predictors Divergence
        Get_PredictorsDivergence(MESH);

		// Poisson System Resolution
        Get_GaussSeidel(P1);

		// New Velocities Calculation
		Get_Velocities(MESH, P1);

		if(Step%10 == 0){

			// Communication to global matrix
			P1.SendMatrixToZeroMP(P.Pres, Global.P);
			if (Problema == 2){ P1.SendMatrixToZeroMP(T.Pres, Global.T); } 		
			P1.SendMatrixToZeroMU(U.Fut, Global.U);
			P1.SendMatrixToZeroMV(V.Fut, Global.V);
			P1.SendMatrixToZeroMW(W.Fut, Global.W);

			// Checking Convergence Criteria
			Get_Stop();

			// Print of .VTK files
			if(Rango == 0){

				sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
				POST1.VTK_GlobalScalar3D("DrivenCavity/", "Presion", FileName_1, MESH, Global.P);
		
				if (Problema == 2){
					POST1.Get_NusseltResults(MESH, Global.T);

					sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
					POST1.VTK_GlobalScalar3D("DrivenCavity/", "Temperatura", FileName_1, MESH, Global.T);
				}
				
				sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
				POST1.VTK_GlobalVectorial3D("DrivenCavity/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);

				// Current Simulation Status
				cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiffGlobal<<endl;	
			}
		}

		// Fields Update
		Get_Update();
		
	}
	
	// Solver Completed
	if(Rango == 0){
		if (Problema == 2){ POST1.Get_NusseltResults(MESH, Global.T); }
		POST1.Get_VelocityResults(MESH, Global.U, Global.V);
		cout<<"Solver Completed."<<endl;
	}

	// Memory Delete

		// ReadData Memory
		R1.Delete_ReadDataMemory();

		// Parallel Memory
		P1.Delete_ParallelMemory();

		// Mesher Memory
		MESH.Delete_MesherMemory();

		// PostProcessing Memory
		POST1.Delete_PostProcessingMemory();

		// Solver
		Delete_VelocityMemory(U);
		Delete_VelocityMemory(V);
		Delete_VelocityMemory(W);

		Delete_EnergyMemory(T);

		Delete_PoissonMemory(A);

		Delete_SolverMemory();
		
	if (Rango == 0){
		cout<<"Memory Deleted."<<endl;
	}

}
