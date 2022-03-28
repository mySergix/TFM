//------------------------------------------------------------------------------------------------//
//                       CPP FILE FOR SPECIES MEMORY FUNCTIONS                                    //
//------------------------------------------------------------------------------------------------//

// Function to allocate memory for each specie on the solver
void Solver::Allocate_StructSpecies(Memory M1){
int n;

	// Species k = 1 to k = N - 1
    for (n = 0; n < N_Species - 1; n++){
		
		// Mass Fraction fields maps in each time step
    	Species[n].Y_Pres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*HP, NZ + 2*HP, 1);
    	Species[n].Y_Fut = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*HP, NZ + 2*HP, 1);

		Species[n].X = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*HP, NZ + 2*HP, 1);

    	// Equation terms
    	Species[n].Convective = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * HP, NY + 2*HP, NZ + 2*HP, 1);
    	Species[n].Diffusive = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * HP, NY + 2*HP, NZ + 2*HP, 1);

		Species[n].ContributionPast = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * HP, NY + 2*HP, NZ + 2*HP, 1);
		Species[n].ContributionPres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * HP, NY + 2*HP, NZ + 2*HP, 1);

		Species[n].D_AlphaMix = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2 * HP, NY + 2*HP, NZ + 2*HP, 1);

    	if (Rango == 0){
			Species[n].Left = M1.AllocateDouble(1, NY, NZ, 1);
    	}
    	else if (Rango == Procesos - 1){
        	Species[n].Right = M1.AllocateDouble(1, NY, NZ, 1);
    	}

    	Species[n].Bottom = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, 1, NZ, 1);
    	Species[n].Top = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, 1, NZ, 1);

    	Species[n].Here = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, NY, 1, 1);
    	Species[n].There = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2, NY, 1, 1);
          
		// JANAF Terms
    	Species[n].Cp_coeff = M1.AllocateDouble(10, 1, 1, 1);
    	Species[n].h_coeff = M1.AllocateDouble(12, 1, 1, 1);
    	Species[n].mu_coeff = M1.AllocateDouble(8, 1, 1, 1);
    	Species[n].lambda_coeff = M1.AllocateDouble(8, 1, 1, 1);

		// Chemical Terms
		Species[n].wk  = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*HP, NZ + 2*HP, 1);

    	if (Rango == 0){
      		Species[n].Global = M1.AllocateDouble(NX + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);
    	}
	
	}

	// Species k = N

		// Mass Fraction fields maps in each time step
    	Species[N_Species-1].Y_Pres = M1.AllocateDouble(Fx[Rango] - Ix[Rango] + 2*HP, NY + 2*HP, NZ + 2*HP, 1);

		if (Rango == 0){
      		Species[N_Species-1].Global = M1.AllocateDouble(NX + 2*HP, NY + 2*Halo, NZ + 2*Halo, 1);
    	}

		// JANAF Terms
    	Species[N_Species-1].Cp_coeff = M1.AllocateDouble(10, 1, 1, 1);
    	Species[N_Species-1].h_coeff = M1.AllocateDouble(12, 1, 1, 1);
    	Species[N_Species-1].mu_coeff = M1.AllocateDouble(8, 1, 1, 1);
    	Species[N_Species-1].lambda_coeff = M1.AllocateDouble(8, 1, 1, 1);
		
}

// Function to delete all the memory allocated for species transport
void Solver::Delete_StructSpecies(){
int n;

	// Species k = 1 to k = N - 1
    for (n = 0; n < N_Species - 1; n++){

		// Mass Fraction fields maps in each time step
    	delete[] Species[n].Y_Pres;
    	delete[] Species[n].Y_Fut;

    	// Equation terms
    	delete[] Species[n].Convective;
    	delete[] Species[n].Diffusive;

		delete[] Species[n].D_ab;
    	
    	if (Rango == 0){
			delete[] Species[n].Left;
    	}
    	else if (Rango == Procesos - 1){
        	delete[] Species[n].Right;
    	}

    	delete[] Species[n].Bottom;
    	delete[] Species[n].Top;

    	delete[] Species[n].Here;
    	delete[] Species[n].There;
          
    	if (Rango == 0){
      		delete[] Species[n].Global;
    	}

	}

	// Species k = N

		// Mass Fraction fields maps in each time step
    	delete[] Species[N_Species-1].Y_Pres;

		if (Rango == 0){
      		delete[] Species[N_Species-1].Global;
    	}

}