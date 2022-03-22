//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR SOLVER EXECUTION                                      //
//------------------------------------------------------------------------------------------------//

// Function to run the solver of the simulation
void CFD_Solver::RunSolver(Memory M1, Parallel P1, Mesher MESH, PostProcessing POST1){	
int i, j, k;

	double Time = 0.0;
	int Step = 0;

	char FileName_1[300]; 

	MaxDiffGlobal = 2.0*ConvergenciaGlobal;

	// Memory Allocation for Navier Stokes Equations
	Allocate_VelocitiesMemory(M1);
    Allocate_PressureMemory(M1);
    Allocate_PoissonCoeffsMemory(M1);
	Allocate_EnergyMemory(M1);
	if (Rango == 0){ Allocate_GlobalMemory(M1); }

    // Initial settings and calculations
    Get_InitialConditions(MESH);
    Get_PoissonCoefficients(MESH);	
    Get_InitialBoundaryConditions(MESH);
	Get_StaticHalos();

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
				
		sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("DrivenCavity/", "Temperatura", FileName_1, MESH, Global.T);

		sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
		POST1.VTK_GlobalVectorial3D("DrivenCavity/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);
	}

	auto InicioLoop = std::chrono::high_resolution_clock::now();
	double TimeRHS = 0.0;
	double TimeKSP = 0.0;
	while(Step < 500){
		//MaxDiffGlobal >= ConvergenciaGlobal
		// New Step
		Step++;
		
		// Update boundary conditions
		Get_UpdateBoundaryConditions(MESH);
        Get_UpdateHalos();
		
		// Halos communication
		P1.CommunicateDataLU(U.Pres, U.Pres);
		P1.CommunicateDataLV(V.Pres, V.Pres);
		P1.CommunicateDataLW(W.Pres, W.Pres);
		
		// Step Time Calculation
		Get_StepTime(MESH, P1);
		Time += DeltaT;
        
		// Diffusion Terms Calculation
        Get_DiffusionU(MESH);
        Get_DiffusionV(MESH);
        Get_DiffusionW(MESH);
        
		// Convective Terms Calculation
		Get_ConvectionU(MESH);
		Get_ConvectionV(MESH);
		Get_ConvectionW(MESH);

		Get_BoussinesqV(MESH);
		
		// Predictors Velocities and Divergence
        Get_ContributionsPredictors();
        Get_PredictorsDivergence(MESH);

		// Poisson System Resolution
       

		auto InicioKSP = std::chrono::high_resolution_clock::now();
			
			 Get_GaussSeidel(P1);

			auto FinalKSP = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> T_KSP = FinalKSP - InicioKSP;
			TimeKSP = TimeKSP + T_KSP.count();

		// New Velocities Calculation
		Get_Velocities(MESH, P1);

		// Energy Equation Calculation

		// Temperature Halo Communication
        P1.CommunicateDataLP(T.Pres, T.Pres);

		Get_DiffusionEnergy(MESH);
		Get_ConvectionEnergy(MESH);
		Get_Temperature();

		if(Step%100 == 0){

			// Communication to global matrix
			P1.SendMatrixToZeroMP(P.Pres, Global.P);
			P1.SendMatrixToZeroMP(T.Pres, Global.T);
			P1.SendMatrixToZeroMU(U.Fut, Global.U);
			P1.SendMatrixToZeroMV(V.Fut, Global.V);
			P1.SendMatrixToZeroMW(W.Fut, Global.W);

			// Checking Convergence Criteria
			Get_Stop();

			// Print of .VTK files
			if(Rango == 0){

				POST1.Get_NusseltResults(MESH, Global.T);

				sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
				POST1.VTK_GlobalScalar3D("DrivenCavity/", "Presion", FileName_1, MESH, Global.P);
		
				sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
				POST1.VTK_GlobalScalar3D("DrivenCavity/", "Temperatura", FileName_1, MESH, Global.T);
				
				sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
				POST1.VTK_GlobalVectorial3D("DrivenCavity/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);

				// Current Simulation Status
				cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiffGlobal<<endl;	
			}
		}

		// Fields Update
		Get_Update();
		
	}
	auto FinalLoop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = FinalLoop - InicioLoop;
	//std::cout << "Elapsed time: " << elapsed.count() << " s\n";

	if (Rango == 0){
		cout<<"Tiempos: "<<endl;
		cout<<"Vector RHS: "<<TimeRHS<<endl;
		cout<<"Solver KSP: "<<TimeKSP<<endl;
		cout<<"Tiempo total: "<<elapsed.count()<<endl;
	}
	// Solver Completed
	if(Rango == 0){
		POST1.Get_NusseltResults(MESH, Global.T);
		POST1.Get_VelocityResults(MESH, Global.U, Global.V);
		cout<<"Solver Completed"<<endl;
	}

}
