//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR SOLVER EXECUTION                                      //
//------------------------------------------------------------------------------------------------//

// Function to run the solver of the simulation
void Solver::RunSolver(Memory M1, Parallel P1, Mesher MESH, PostProcessing POST1){	
int i, j, k;

	PetscMPIInt rank, size;
	PetscInt	m = (Fx[Rango] - Ix[Rango]) * (NY) * (NZ); // Number of local rows
	PetscInt	n = (Fx[Rango] - Ix[Rango]) * (NY) * (NZ); // Number of local columns

	PetscInt	M = (NX) * (NY) * (NZ); // Number of global rows
	PetscInt	N = (NX) * (NY) * (NZ); // Number of global columns

	KSP         ksp;
 	PC          pc;
 	PetscBool   pinPressure;

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

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
	NNZ = 0;
	Get_NonZero_NumberElements();

	PetscMalloc1(NNZ, &Val_Laplacian);  
    PetscMalloc1(NNZ, &Col_Ind);
    PetscMalloc1((Fx[Rango] - Ix[Rango]) * (NY) * (NZ) + 1, &Row_Ptr);
	PetscMalloc1((Fx[Rango] - Ix[Rango]) * (NY) * (NZ), &RHS_Ind);  
	PetscMalloc1((Fx[Rango] - Ix[Rango]) * (NY) * (NZ), &RHS);  
/*
	VecCreate(PETSC_COMM_WORLD, &B_RHS);
	VecSetType(B_RHS, VECMPI);
	VecCreate(PETSC_COMM_WORLD, &X_Sol);
	VecSetType(X_Sol, VECMPI);

	VecSetSizes(B_RHS, n, M);
	VecSetSizes(X_Sol, n, N);
	*/
	Get_CSR_LaplacianMatrix();
	Get_RHS_VectorIndex();

	MatCreate(PETSC_COMM_WORLD, &A_Matrix);

	MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, m, n, PETSC_DECIDE, PETSC_DECIDE, Row_Ptr, Col_Ind, Val_Laplacian, &A_Matrix);
	MatAssemblyBegin(A_Matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A_Matrix, MAT_FINAL_ASSEMBLY);

	PetscInt m_matriz, n_matriz;

	MatGetSize(A_Matrix, &m_matriz, &n_matriz);

	if (Rango == 1){
		cout<<"m_matriz: "<<m_matriz<<endl;
		cout<<"n_matriz: "<<n_matriz<<endl;
	}

	// Using PETSc
	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetOperators(ksp, A_Matrix, A_Matrix);


	//KSPSetType(ksp, KSPCG); // 19 s aprox
	//KSPSetType(ksp, KSPGROPPCG); // 19 s. aprox
	//KSPSetType(ksp, KSPPIPECG); // 20 s aprox
	//KSPSetType(ksp, KSPPIPEPRCG); // 30 s aprox
	//KSPSetType(ksp, KSPPIPECG2); // 22 s aprox
	KSPSetType(ksp, KSPCG);

	//KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
	KSPGetPC(ksp, &pc);
	//KSPSetPCSide(ksp, PC_LEFT);
	//PCSetType(pc, PCJACOBI);
	KSPSetTolerances(ksp, 1e-2, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPSetFromOptions(ksp);
	KSPSetUp(ksp);
	//KSPGetPC(ksp, &pc);
	//PCSetType(pc, PCJACOBI);
	//KSPSetTolerances(ksp, 1.e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

	

	VecCreate(PETSC_COMM_WORLD, &X_Sol);
	VecSetSizes(X_Sol, m, M);
	VecSetFromOptions(X_Sol);
	PetscObjectSetName((PetscObject)X_Sol,"Approx. Solution");
	VecDuplicate(X_Sol, &B_RHS);
	PetscObjectSetName((PetscObject)B_RHS,"Right hand side");
	//VecSet(X_Sol, 0.0);
	//VecSet(B_RHS, 0.0);
	//MatMult(A_Matrix, X_Sol, B_RHS);

	//MatView(A_Matrix, 	PETSC_VIEWER_STDOUT_SELF 	);
/*
	VecCreate(PETSC_COMM_WORLD, &B_RHS);
 	VecSetSizes(B_RHS, m, PETSC_DECIDE);
	VecSetFromOptions(B_RHS);
	VecDuplicate(B_RHS, &X_Sol);

	VecSet(X_Sol,1.0);
    MatMult(A_Matrix, X_Sol, B_RHS);

	*/

	//
	
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
	//

	auto InicioLoop = std::chrono::high_resolution_clock::now();
	double TimeRHS = 0.0;
	double TimeKSP = 0.0;
	while(MaxDiffGlobal >= ConvergenciaGlobal){
		//
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

		if (Problema == 2){	Get_BoussinesqV(MESH); }
			
		// Predictors Velocities and Divergence
        Get_ContributionsPredictors();
        Get_PredictorsDivergence(MESH);

		// Poisson System Resolution
        //Get_GaussSeidel(P1);

			auto InicioRHS = std::chrono::high_resolution_clock::now();
			// Using PETSc
			VecSetValues(B_RHS, (Fx[Rango] - Ix[Rango]) * (NY) * (NZ), (const PetscInt*)RHS_Ind, (const PetscScalar*)RHS, INSERT_VALUES);
			VecAssemblyBegin(B_RHS);
			VecAssemblyEnd(B_RHS);
			
			auto FinalRHS = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> T_RHS = FinalRHS - InicioRHS;
			TimeRHS = TimeRHS + T_RHS.count();

			auto InicioKSP = std::chrono::high_resolution_clock::now();
			//VecView(B_RHS, 	PETSC_VIEWER_STDOUT_WORLD);
			KSPSolve(ksp, B_RHS, X_Sol);

			auto FinalKSP = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> T_KSP = FinalKSP - InicioKSP;
			TimeKSP = TimeKSP + T_KSP.count();

			VecGetArray(X_Sol, &X_Sol_Array);
			Get_LocalSolution();

		// New Velocities Calculation
		Get_Velocities(MESH, P1);

		// Energy Equation Calculation
		if (Problema == 2){
			// Temperature Halo Communication
        	P1.CommunicateDataLP(T.Pres, T.Pres);

			Get_DiffusionEnergy(MESH);
			Get_ConvectionEnergy(MESH);
			Get_Temperature();
		}
		
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
		//POST1.Get_NusseltResults(MESH, Global.T);
		POST1.Get_VelocityResults(MESH, Global.U, Global.V);
		cout<<"Solver Completed"<<endl;
	}

	// Petsc Memory free
	PetscFree(Val_Laplacian);
	PetscFree(Col_Ind);
	PetscFree(Row_Ptr);

	KSPDestroy(&ksp);
	VecDestroy(&B_RHS);
	VecDestroy(&X_Sol);

	MatDestroy(&A_Matrix);
}
