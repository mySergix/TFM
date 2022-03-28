//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR SOLVER EXECUTION                                      //
//------------------------------------------------------------------------------------------------//

// Function to run the solver of the simulation
void Solver::RunSolver(Memory M1, ReadData R1, Parallel P1, Mesher MESH, PostProcessing POST1){	
int i, j, k;

	// PETSc Library
	PetscMPIInt rank, size;
	PetscInt	m = (Fx[Rango] - Ix[Rango]) * (NY) * (NZ); // Number of local rows
	PetscInt	n = (Fx[Rango] - Ix[Rango]) * (NY) * (NZ); // Number of local columns

	PetscInt	M = (NX) * (NY) * (NZ); // Number of global rows
	PetscInt	N = (NX) * (NY) * (NZ); // Number of global columns
 
	KSP         ksp; // Krylov Subspace Solver
 	PC          pc; // Preconditioner
	MatNullSpace nullspace_Amatrix; // Nullspace setting for Singular Laplacian Matrix

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	double Time = 0.0;
	int Step = 0;

	char FileName_1[300]; 

	MaxDiffGlobal = 2.0*ConvergenciaGlobal;

	// Runge Kutta 3rd Order Coefficients
	Get_RK_Coefficients(RK);

	// Memory Allocation for Navier Stokes Equations
	Allocate_VelocitiesMemory(M1);
    Allocate_PressureMemory(M1);
    Allocate_PoissonCoeffsMemory(M1);
	Allocate_EnergyMemory(M1);
	Allocate_StructSpecies(M1);
	if (Rango == 0){ Allocate_GlobalMemory(M1); }

	// Species Data Reading
	Read_SpeciesName("Species_Data.txt");
	Read_AllSpeciesData();

	/*
	// Laplacian Matrix Calculations
	Get_PoissonCoefficients(MESH);	
	NNZ = 0; // Non-Zero Elements
	Get_NonZero_NumberElements();

	// Memory Allocation for CSR Laplacian Matrix and Linear System
	PetscMalloc1(NNZ, &Val_Laplacian);  
    PetscMalloc1(NNZ, &Col_Ind);
	PetscMalloc1(m + 1, &Row_Ptr);
	PetscMalloc1(m, &RHS_Ind);  
	PetscMalloc1(m, &RHS);  

	Get_CSR_LaplacianMatrix();
	Get_RHS_VectorIndex();

	MatCreate(PETSC_COMM_WORLD, &A_Matrix);

	MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, m, n, PETSC_DECIDE, PETSC_DECIDE, Row_Ptr, Col_Ind, Val_Laplacian, &A_Matrix);
	MatAssemblyBegin(A_Matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A_Matrix, MAT_FINAL_ASSEMBLY);

	// Krylov Subspace Solver Creation
	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetOperators(ksp, A_Matrix, A_Matrix);
	KSPSetType(ksp, KSPCG);
	KSPGetPC(ksp, &pc);
	KSPSetTolerances(ksp, 1e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
	MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace_Amatrix);
	MatSetNullSpace(A_Matrix, nullspace_Amatrix);
	MatSetTransposeNullSpace(A_Matrix, nullspace_Amatrix);
	KSPSetFromOptions(ksp);
	KSPSetUp(ksp);

	// RHS and Solution Vector Creation
	VecCreate(PETSC_COMM_WORLD, &X_Sol); 
	VecSetSizes(X_Sol, m, M);
	VecSetFromOptions(X_Sol);
	PetscObjectSetName((PetscObject)X_Sol, "Pressure Solution");
	VecDuplicate(X_Sol, &B_RHS);
	PetscObjectSetName((PetscObject)B_RHS,"Right Hand Side");

    // Initial Simulation settings
	Get_InitialConditions(MESH);
    Get_StaticBoundaryConditions_Velocities(MESH);
	Get_StaticHalos_Velocity(U.Pres, V.Pres, W.Pres);
	//Get_MassFlowValue(MESH);

	Get_UpdateBoundaryConditions_Velocities(U.Pres, V.Pres, W.Pres);
	Get_UpdateHalos_Velocity(U.Pres, V.Pres, W.Pres);

	Get_StaticBoundaryConditions_Temperatures(MESH);
	Get_StaticHalos_Temperatures(T.Pres);

	Get_UpdateBoundaryConditions_Temperatures(MESH, T.Pres);
	Get_UpdateHalos_Temperatures(T.Pres);

	Get_Species_InitialConditions();
	Get_Species_StaticBoundaryConditions(MESH);
	Get_Species_StaticHalos();

	Get_Species_UpdateBoundaryConditions(MESH);
	Get_Species_UpdateHalos();
	Get_Species_DiffusionCoefficients();
	Get_DiffusiveTimeStep(MESH);

	// Communication to global matrix at step 0
	P1.SendMatrixToZeroMP(P.Pres, Global.P);
	P1.SendMatrixToZeroMU(U.Fut, Global.U);
	P1.SendMatrixToZeroMV(V.Fut, Global.V);
	P1.SendMatrixToZeroMW(W.Fut, Global.W);
	P1.SendMatrixToZeroMP(T.Pres, Global.T);
	P1.SendMatrixToZeroMP(Species[0].Y_Pres, Species[0].Global);

	// Print of Step 0 .VTK files
	if(Rango == 0){
		sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("ChannelTransport/", "Presion", FileName_1, MESH, Global.P);

		sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("ChannelTransport/", "Temperatura", FileName_1, MESH, Global.T);
		
		sprintf(FileName_1, "MapaFracciones_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("ChannelTransport/", "Fracciones", FileName_1, MESH, Species[0].Global);

		sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
		POST1.VTK_GlobalVectorial3D("ChannelTransport/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);
	}

	// Time Variables
	double Time_RK = 0.0;
	double Time_PoissonSolver = 0.0;

	// Initial Loop Time
	auto InicioLoop = std::chrono::high_resolution_clock::now();

	while(MaxDiffGlobal >= ConvergenciaGlobal){
		
		// New Step
		Step++;

		// Step Time Calculation
		Get_StepTime(MESH);
		Time += DeltaT;
        
		// RK3 Integration (Velocity + Temperature)
		auto Inicio_RK = std::chrono::high_resolution_clock::now();
		Get_RK_Integration(MESH, P1);
		auto Final_RK = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> T_RK = Final_RK - Inicio_RK;
		Time_RK += T_RK.count();

		// Predictors Divergence
        Get_PredictorsDivergence(MESH);

		// RHS Vector Assembly
		VecSetValues(B_RHS, m, (const PetscInt*)RHS_Ind, (const PetscScalar*)RHS, INSERT_VALUES);
		VecAssemblyBegin(B_RHS);
		VecAssemblyEnd(B_RHS);
			
		// Poisson System Resolution
        auto Inicio_PoissonSolver = std::chrono::high_resolution_clock::now();
		KSPSolve(ksp, B_RHS, X_Sol);
		auto Final_PoissonSolver = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> T_PoissonSolver = Final_PoissonSolver - Inicio_PoissonSolver;
		Time_PoissonSolver += T_PoissonSolver.count();

		VecGetArray(X_Sol, &X_Sol_Array);
		Get_LocalSolution();

		// New Velocities Calculation
		Get_Velocities(MESH, P1);

		// New Temperatures Calculation
		Get_UpdateBoundaryConditions_Temperatures(MESH, T.Pres);
		Get_UpdateHalos_Temperatures(T.Pres);

		P1.CommunicateDataLP(T.Pres, T.Pres);
		Get_DiffusionEnergy(MESH, T.Pres);
		Get_ConvectionEnergy(MESH, T.Pres, U.Pres, V.Pres, W.Pres);
		Get_EnergyContributions();
        Get_NewTemperatures(MESH);

		// New Species Mass Fraction Calculation
		Get_Species_UpdateBoundaryConditions(MESH);
		Get_Species_UpdateHalos();

		P1.CommunicateDataLP(Species[0].Y_Pres, Species[0].Y_Pres);
		Get_Species_Diffusion(MESH);
		Get_Species_Convection(MESH);

		Get_SpeciesContributions();
		Get_Species_MassFraction();
		Get_Species_MassConservation();

		

		if (Step%100 == 0){
			// Checking Convergence Criteria
			Get_Stop();

			if (Rango == 0){
				// Current Simulation Status
				cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiffGlobal<<", DeltaT: "<<DeltaT<<endl;	
			}
		}

		// Fields Update
		Get_Update();

		if(Step%200 == 0){

			// Communication to global matrix
			P1.SendMatrixToZeroMP(P.Pres, Global.P);
			P1.SendMatrixToZeroMU(U.Pres, Global.U);
			P1.SendMatrixToZeroMV(V.Pres, Global.V);
			P1.SendMatrixToZeroMW(W.Pres, Global.W);
			P1.SendMatrixToZeroMP(T.Pres, Global.T);
			P1.SendMatrixToZeroMP(Species[0].Y_Pres, Species[0].Global);
			P1.SendBoCoToZeroMP(T.Bottom, Global.Tbottom);
			P1.SendBoCoToZeroMP(Species[0].Bottom, Global.Ybottom);

			// Print of .VTK files
			if(Rango == 0){
				
				POST1.Get_NusseltNumber(MESH, Global.T, Global.U, Global.Tbottom, Twater, Dh, Qs);
				POST1.Get_SherwoodNumber(MESH, Species[0].Global, Global.U, Global.Ybottom, Dh);

				sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
				POST1.VTK_GlobalScalar3D("ChannelTransport/", "Presion", FileName_1, MESH, Global.P);
		
				sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
				POST1.VTK_GlobalScalar3D("ChannelTransport/", "Temperatura", FileName_1, MESH, Global.T);
				
				sprintf(FileName_1, "MapaFracciones_Step_%d", Step);
				POST1.VTK_GlobalScalar3D("ChannelTransport/", "Fracciones", FileName_1, MESH, Species[0].Global);

				sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
				POST1.VTK_GlobalVectorial3D("ChannelTransport/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);
		
			}
		}
	}

	// Final Loop Time
	auto FinalLoop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> Time_Loop = FinalLoop - InicioLoop;

	// Communication to global matrix
	P1.SendMatrixToZeroMP(P.Pres, Global.P);
	P1.SendMatrixToZeroMU(U.Fut, Global.U);
	P1.SendMatrixToZeroMV(V.Fut, Global.V);
	P1.SendMatrixToZeroMW(W.Fut, Global.W);
	P1.SendMatrixToZeroMP(T.Pres, Global.T);
	P1.SendMatrixToZeroMP(Species[0].Y_Pres, Species[0].Global);
	P1.SendBoCoToZeroMP(T.Bottom, Global.Tbottom);
	P1.SendBoCoToZeroMP(Species[0].Bottom, Global.Ybottom);

	// Print of .VTK files
	if(Rango == 0){	

		sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("ChannelTransport/", "Presion", FileName_1, MESH, Global.P);
		
		sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("ChannelTransport/", "Temperatura", FileName_1, MESH, Global.T);
				
		sprintf(FileName_1, "MapaFracciones_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("ChannelTransport/", "Fracciones", FileName_1, MESH, Species[0].Global);

		sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
		POST1.VTK_GlobalVectorial3D("ChannelTransport/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);

		// Post Processing Calculations
		POST1.Get_NusseltNumber(MESH, Global.T, Global.U, Global.Tbottom, Twater, Dh, Qs);
		POST1.Get_SherwoodNumber(MESH, Species[0].Global, Global.U, Global.Ybottom, Dh);

		// Solver Completed
		cout<<"Solver Completed"<<endl;
		cout<<endl;
		cout<<"Solver Times:"<<endl;
		cout<<"RK3 Integration: "<<Time_RK<<" s, "<<100 * (Time_RK / Time_Loop.count())<<" %"<<endl;
		cout<<"Poisson Solver: "<<Time_PoissonSolver<<" s, "<<100 * (Time_PoissonSolver / Time_Loop.count())<<" %"<<endl;
		cout<<"Loop Time: "<<Time_Loop.count()<<" s"<<endl;
	}

	// Petsc Memory free
	PetscFree(Val_Laplacian);
	PetscFree(Col_Ind);
	PetscFree(Row_Ptr);

	PetscFree(RHS_Ind);
	PetscFree(RHS);

	KSPDestroy(&ksp);
	VecDestroy(&B_RHS);
	VecDestroy(&X_Sol);

	MatDestroy(&A_Matrix);

	// Classes Memory Delete

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
		Delete_StructSpecies();

		Delete_SolverMemory();
	*/	
	if (Rango == 0){
		cout<<endl;
		cout<<"Memory Deleted."<<endl;
	}

}