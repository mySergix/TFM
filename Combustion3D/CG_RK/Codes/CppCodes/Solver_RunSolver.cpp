//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR SOLVER EXECUTION                                      //
//------------------------------------------------------------------------------------------------//

// Function to run the solver of the simulation
void Solver::RunSolver(Memory M1, ReadData R1, Parallel P1, Mesher MESH, PostProcessing POST1){	
int i, j, k;

	// PETSc Library
	PetscMPIInt rank, size;
	PetscInt	m = ProcessNodesP; // Number of local rows
	PetscInt	n = ProcessNodesP; // Number of local columns

	PetscInt	M = TotalNodesP; // Number of global rows
	PetscInt	N = TotalNodesP; // Number of global columns
 
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

	// Processes Internal Walls Check
	Get_ProcessesInternalCheck();

	// Memory Allocation for Navier Stokes Equations
	Allocate_VelocitiesMemory(M1);
    Allocate_PressureMemory(M1);
    Allocate_PoissonCoeffsMemory(M1);
	Allocate_EnergyMemory(M1);
	if (Rango == 0){ Allocate_GlobalMemory(M1); }
	
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
	
	Get_CSR_LaplacianMatrix(MESH);
	Get_RHS_VectorIndex(MESH);
	
	MatCreate(PETSC_COMM_WORLD, &A_Matrix);
	
	MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD, m, n, PETSC_DECIDE, PETSC_DECIDE, Row_Ptr, Col_Ind, Val_Laplacian, &A_Matrix);
	MatAssemblyBegin(A_Matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A_Matrix, MAT_FINAL_ASSEMBLY);
	
	// Krylov Subspace Solver Creation
	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetOperators(ksp, A_Matrix, A_Matrix);
	KSPSetType(ksp, KSPCG);
	KSPGetPC(ksp, &pc);
	KSPSetTolerances(ksp, 1e-4, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
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

	Get_StaticHalos_Velocity(MESH, U.Pres, V.Pres, W.Pres);
	Get_StaticHalos_Velocity(MESH, U.Fut, V.Fut, W.Fut);
	Get_StaticHalos_Velocity(MESH, U.New_Velocity, V.New_Velocity, W.New_Velocity);
	
	Get_UpdateBoundaryConditions_Velocities(MESH, U.Pres, V.Pres, W.Pres);
	Get_UpdateBoundaryConditions_Velocities(MESH, U.Fut, V.Fut, W.Fut);

	Get_UpdateHalos_Velocity(MESH, U.Pres, V.Pres, W.Pres);
	Get_UpdateHalos_Velocity(MESH, U.Fut, V.Fut, W.Fut);

	Get_StaticBoundaryConditions_Temperatures(MESH);
	Get_StaticHalos_Temperatures(MESH, T.Pres);
	
	Get_DiffusiveTimeStep(MESH);

	// Communication to global matrix at step 0
	P1.SendMatrixToZeroMP(P.Pres, Global.P);
	P1.SendMatrixToZeroMU(U.Fut, Global.U);
	P1.SendMatrixToZeroMV(V.Fut, Global.V);
	P1.SendMatrixToZeroMW(W.Fut, Global.W);
	P1.SendMatrixToZeroMP(T.Pres, Global.T); 

	// Print of Step 0 .VTK files
	if(Rango == 0){
		sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("DrivenCavity/", "Presion", FileName_1, MESH, Global.P);

		sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("DrivenCavity/", "Temperatura", FileName_1, MESH, Global.T);
		
		sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
		POST1.VTK_GlobalVectorial3D("DrivenCavity/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);
	}
	
	while(MaxDiffGlobal >= ConvergenciaGlobal){
		
		// New Step
		Step++;

		// Step Time Calculation
		Get_StepTime(MESH);
		Time += DeltaT;
        
		// RK Integration (Velocity Predictor)
		Get_RK_Integration(MESH, P1);
		
		// Predictors Divergence
        Get_PredictorsDivergence(MESH);
		
		// RHS Vector Assembly
		VecSetValues(B_RHS, m, (const PetscInt*)RHS_Ind, (const PetscScalar*)RHS, INSERT_VALUES);
		VecAssemblyBegin(B_RHS);
		VecAssemblyEnd(B_RHS);
		
		// Poisson System Resolution
        
		KSPSolve(ksp, B_RHS, X_Sol);

		VecGetArray(X_Sol, &X_Sol_Array);
		Get_LocalSolution(MESH);
		
		// New Velocities Calculation
		Get_Velocities(MESH, P1);
		
		/*
		// New Temperatures Calculation
		Get_UpdateBoundaryConditions_Temperatures(MESH, T.Pres);
		Get_UpdateHalos_Temperatures(MESH, T.Pres);

		P1.CommunicateDataLP(T.Pres, T.Pres);
		Get_DiffusionEnergy(MESH, T.Pres);
		Get_ConvectionEnergy(MESH, T.Pres, U.Pres, V.Pres, W.Pres);
		Get_EnergyContributions(MESH);
        Get_NewTemperatures(MESH);
		*/

		if (Step%10 == 0){
			// Checking Convergence Criteria
			Get_Stop(MESH);

			if (Rango == 0){

				// Current Simulation Status
				cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiffGlobal<<", DeltaT: "<<DeltaT<<endl;	
			}

		}
		
		// Fields Update
		Get_Update(MESH);
		
		if(Step%200 == 0){
			
			// Pressure Halos
			Get_PressureHalos();

			// Communication to global matrix
			P1.SendMatrixToZeroMP(P.Pres, Global.P);
			P1.SendMatrixToZeroMU(U.Fut, Global.U);
			P1.SendMatrixToZeroMV(V.Fut, Global.V);
			P1.SendMatrixToZeroMW(W.Fut, Global.W);
			P1.SendMatrixToZeroMP(T.Pres, Global.T);

			// Print of .VTK files
			if(Rango == 0){
				
				POST1.Get_GlobalScalarHalos(MESH, Global.P);
				POST1.Get_GlobalVectorialHalos(Global.U, Global.V, Global.W);

				sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
				POST1.VTK_GlobalScalar3D("DrivenCavity/", "Presion", FileName_1, MESH, Global.P);
		
				sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
				POST1.VTK_GlobalScalar3D("DrivenCavity/", "Temperatura", FileName_1, MESH, Global.T);
				
				sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
				POST1.VTK_GlobalVectorial3D("DrivenCavity/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);
		
			}
			
		}
		
	}
	

	// Communication to global matrix
	P1.SendMatrixToZeroMP(P.Pres, Global.P);
	P1.SendMatrixToZeroMU(U.Fut, Global.U);
	P1.SendMatrixToZeroMV(V.Fut, Global.V);
	P1.SendMatrixToZeroMW(W.Fut, Global.W);
	P1.SendMatrixToZeroMP(T.Pres, Global.T);

	// Print of .VTK files
	if(Rango == 0){		
		
		POST1.Get_GlobalScalarHalos(MESH, Global.P);
		POST1.Get_GlobalVectorialHalos(Global.U, Global.V, Global.W);
		
		sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("DrivenCavity/", "Presion", FileName_1, MESH, Global.P);
		
		sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
		POST1.VTK_GlobalVectorial3D("DrivenCavity/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);

		sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
		POST1.VTK_GlobalScalar3D("DrivenCavity/", "Temperatura", FileName_1, MESH, Global.T);

		// Post Processing Calculations
		for (j = 0; j < NY_1 + NY_2; j++){
			cout<<Global.P[GP(NX_1-2,j,NZ/2,0)]<<", "<<Global.P[GP(NX_1-1,j,NZ/2,0)]<<endl;
		}

		// Solver Completed
		cout<<"Solver Completed"<<endl;
		cout<<endl;
		
	}
	
/*
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
	*/
/*
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
		Delete_SolverMemory();
	*/	
	if (Rango == 0){
		cout<<endl;
		cout<<"Memory Deleted."<<endl;
	}

}