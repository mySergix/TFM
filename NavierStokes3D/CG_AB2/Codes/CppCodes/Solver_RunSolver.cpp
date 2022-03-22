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

	if (Problema == 2){
		Get_StaticBoundaryConditions_Temperatures(MESH);
		Get_StaticHalos_Temperatures(T.Pres);
	}

	Get_DiffusiveTimeStep(MESH);

	// Communication to global matrix at step 0
	P1.SendMatrixToZeroMP(P.Pres, Global.P);
	P1.SendMatrixToZeroMU(U.Fut, Global.U);
	P1.SendMatrixToZeroMV(V.Fut, Global.V);
	P1.SendMatrixToZeroMW(W.Fut, Global.W);
	if (Problema == 2){ P1.SendMatrixToZeroMP(T.Pres, Global.T); }

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

	// Time Variables
	double Time_RK3 = 0.0;

	double Time_DiffusionU = 0.0;
	double Time_DiffusionV = 0.0;
	double Time_DiffusionW = 0.0;

	double Time_ConvectionU = 0.0;
	double Time_ConvectionV = 0.0;
	double Time_ConvectionW = 0.0;

	double Time_PoissonSolver = 0.0;
	

	// Initial Loop Time
	auto InicioLoop = std::chrono::high_resolution_clock::now();
	// MaxDiffGlobal >= ConvergenciaGlobal
	while(Step < 5000){
		
		// New Step
		Step++;

		// Step Time Calculation
		Get_StepTime(MESH);
		Time += DeltaT;
        
		// RK3 Integration (Velocity + Temperature)
		auto Inicio_RK3 = std::chrono::high_resolution_clock::now();
		
		// Boundary Conditions Update
        Get_UpdateBoundaryConditions_Velocities(U.Pres, V.Pres, W.Pres);
        Get_UpdateHalos_Velocity(U.Pres, V.Pres, W.Pres);

        // Halos communication
		P1.CommunicateDataLU(U.Pres, U.Pres);
		P1.CommunicateDataLV(V.Pres, V.Pres);
		P1.CommunicateDataLW(W.Pres, W.Pres);
        
        // Diffusion Terms Calculation

		auto Inicio_DiffusionU = std::chrono::high_resolution_clock::now();
        Get_DiffusionU(MESH, U.Pres);
		auto Final_DiffusionU = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> T_DiffusionU = Final_DiffusionU - Inicio_DiffusionU;
		Time_DiffusionU += T_DiffusionU.count();

		auto Inicio_DiffusionV = std::chrono::high_resolution_clock::now();
        Get_DiffusionV(MESH, V.Pres);
		auto Final_DiffusionV = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> T_DiffusionV = Final_DiffusionV - Inicio_DiffusionV;
		Time_DiffusionV += T_DiffusionV.count();

		auto Inicio_DiffusionW = std::chrono::high_resolution_clock::now();
        Get_DiffusionW(MESH, W.Pres);
        auto Final_DiffusionW = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> T_DiffusionW = Final_DiffusionW - Inicio_DiffusionW;
		Time_DiffusionW += T_DiffusionW.count();

		// Convective Terms Calculation

		auto Inicio_ConvectionU = std::chrono::high_resolution_clock::now();
		Get_ConvectionU(MESH, U.Pres, V.Pres, W.Pres);
		auto Final_ConvectionU = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> T_ConvectionU = Final_ConvectionU - Inicio_ConvectionU;
		Time_ConvectionU += T_ConvectionU.count();

		auto Inicio_ConvectionV = std::chrono::high_resolution_clock::now();
		Get_ConvectionV(MESH, U.Pres, V.Pres, W.Pres);
		auto Final_ConvectionV = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> T_ConvectionV = Final_ConvectionV - Inicio_ConvectionV;
		Time_ConvectionV += T_ConvectionV.count();

		auto Inicio_ConvectionW = std::chrono::high_resolution_clock::now();
		Get_ConvectionW(MESH, U.Pres, V.Pres, W.Pres);
		auto Final_ConvectionW = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> T_ConvectionW = Final_ConvectionW - Inicio_ConvectionW;
		Time_ConvectionW += T_ConvectionW.count();

		Get_VelocityContributions();
		Get_PredictorVelocities();

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
		if (Problema == 2){
			Get_UpdateBoundaryConditions_Temperatures(T.Pres);
			Get_UpdateHalos_Temperatures(T.Pres);

			P1.CommunicateDataLP(T.Pres, T.Pres);
			Get_DiffusionEnergy(MESH, T.Pres);
			Get_ConvectionEnergy(MESH, T.Pres, U.Pres, V.Pres, W.Pres);
			Get_EnergyContributions();
            Get_NewTemperatures();
		}

		if (Step%100 == 0){
			// Checking Convergence Criteria
			Get_Stop();

			if (Rango == 0){

				// Current Simulation Status
				cout<<"Step: "<<Step<<", Total time: "<<Time<<", MaxDif: "<<MaxDiffGlobal<<", DeltaT: "<<DeltaT<<endl;	
			}

		}

		if(Step%500 == 0){

			// Communication to global matrix
			P1.SendMatrixToZeroMP(P.Pres, Global.P);
			P1.SendMatrixToZeroMU(U.Fut, Global.U);
			P1.SendMatrixToZeroMV(V.Fut, Global.V);
			P1.SendMatrixToZeroMW(W.Fut, Global.W);
			if (Problema == 2){P1.SendMatrixToZeroMP(T.Pres, Global.T); }

			// Print of .VTK files
			if(Rango == 0){
				
				POST1.Get_VelocityResults(MESH, Global.U, Global.V, K);
				if (Problema == 2){
					POST1.Get_NusseltResults(MESH, Global.T);
				}

				sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
				POST1.VTK_GlobalScalar3D("DrivenCavity/", "Presion", FileName_1, MESH, Global.P);
		
				if (Problema == 2){
					POST1.Get_NusseltResults(MESH, Global.T);

					sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
					POST1.VTK_GlobalScalar3D("DrivenCavity/", "Temperatura", FileName_1, MESH, Global.T);
				}
				
				sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
				POST1.VTK_GlobalVectorial3D("DrivenCavity/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);
		
			}

		}

		// Fields Update
		Get_Update();
		
	}

	// Final Loop Time
	auto FinalLoop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> Time_Loop = FinalLoop - InicioLoop;

	// Communication to global matrix
			P1.SendMatrixToZeroMP(P.Pres, Global.P);
			P1.SendMatrixToZeroMU(U.Fut, Global.U);
			P1.SendMatrixToZeroMV(V.Fut, Global.V);
			P1.SendMatrixToZeroMW(W.Fut, Global.W);
			if (Problema == 2){P1.SendMatrixToZeroMP(T.Pres, Global.T); }

			// Print of .VTK files
			if(Rango == 0){
				
				POST1.Get_VelocityResults(MESH, Global.U, Global.V, K);
				if (Problema == 2){
					POST1.Get_NusseltResults(MESH, Global.T);
				}

				sprintf(FileName_1, "MapaPresiones_Step_%d", Step);
				POST1.VTK_GlobalScalar3D("DrivenCavity/", "Presion", FileName_1, MESH, Global.P);
		
				if (Problema == 2){
					POST1.Get_NusseltResults(MESH, Global.T);

					sprintf(FileName_1, "MapaTemperaturas_Step_%d", Step);
					POST1.VTK_GlobalScalar3D("DrivenCavity/", "Temperatura", FileName_1, MESH, Global.T);
				}
				
				sprintf(FileName_1, "MapaVelocidades_Step_%d", Step);
				POST1.VTK_GlobalVectorial3D("DrivenCavity/", "Velocidades", FileName_1, MESH, Global.U, Global.V, Global.W);
		
			}

	// Solver Completed
	if(Rango == 0){
		POST1.Get_VelocityResults(MESH, Global.U, Global.V, K);
		if (Problema == 2){
			POST1.Get_NusseltResults(MESH, Global.T);
		}
		cout<<"Solver Completed"<<endl;
		cout<<endl;
		cout<<"Solver Times:"<<endl;
		cout<<"DiffusionU: "<<Time_DiffusionU<<" s, "<<100 * (Time_DiffusionU / Time_Loop.count())<<" %"<<endl;
		cout<<"DiffusionV: "<<Time_DiffusionV<<" s, "<<100 * (Time_DiffusionV / Time_Loop.count())<<" %"<<endl;
		cout<<"DiffusionW: "<<Time_DiffusionW<<" s, "<<100 * (Time_DiffusionW / Time_Loop.count())<<" %"<<endl;
		cout<<endl;
		cout<<"ConvectionU: "<<Time_ConvectionU<<" s, "<<100 * (Time_ConvectionU / Time_Loop.count())<<" %"<<endl;
		cout<<"ConvectionV: "<<Time_ConvectionV<<" s, "<<100 * (Time_ConvectionV / Time_Loop.count())<<" %"<<endl;
		cout<<"ConvectionW: "<<Time_ConvectionW<<" s, "<<100 * (Time_ConvectionW / Time_Loop.count())<<" %"<<endl;
		cout<<endl;
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

	// Classes Mmeory Delete

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
		cout<<endl;
		cout<<"Memory Deleted."<<endl;
	}

}