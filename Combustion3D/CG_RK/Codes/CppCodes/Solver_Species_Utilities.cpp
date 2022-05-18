//------------------------------------------------------------------------------------------------//
//                       CPP FILE FOR SPECIES UTILITIES FUNCTIONS                                 //
//------------------------------------------------------------------------------------------------//

// Function to set the initial fields
void Solver::Get_Species_InitialConditions(Mesher MESH){
int i, j, k, n;

    // Species 0 - CH4
    Species[0].InitialY = 0.0;

    // Species 1 - O2
    Species[1].InitialY = 0.21;

    // Species 2 - CO2
    Species[2].InitialY = 0.0;

    // Species 3 - H2O
    Species[3].InitialY = 0.0;

    // Species 4 - N2
    Species[4].InitialY = 0.79;

    for (n = 0; n < N_Species - 1; n++){

        // Fluid Domain
        for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		    for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
			    for(k = 0; k < NZ; k++){	
				    Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].InitialY;
				    Species[n].Y_Fut[LP(i,j,k,0)] = Species[n].InitialY;

                    Species[n].ContributionPres[LP(i,j,k,0)] = 0.0;
                    Species[n].ContributionPast[LP(i,j,k,0)] = 0.0;
			    }
		    }
	    }
        
        // Solid Domain
        for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
            if (MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] > 0){
                for (j = 0; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j++){
                    for(k = 0; k < NZ; k++){	
				        Species[n].Y_Pres[LP(i,j,k,0)] = 0.0;
				        Species[n].Y_Fut[LP(i,j,k,0)] = 0.0;

                        Species[n].ContributionPres[LP(i,j,k,0)] = 0.0;
                        Species[n].ContributionPast[LP(i,j,k,0)] = 0.0;
			        }
                }
            }  
        }

    }
		

    // Species n = N_Species (N2)
    n = N_Species - 1;

    // Fluid Domain
        for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		    for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
			    for(k = 0; k < NZ; k++){	
				    Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].InitialY;
			    }
		    }
	    }
        
        // Solid Domain
        for (i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
            if (MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] > 0){
                for (j = 0; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j++){
                    for(k = 0; k < NZ; k++){	
				        Species[n].Y_Pres[LP(i,j,k,0)] = 0.0;
			        }
                }
            }  
        }

}

// Function to calculate the species diffusion time step
void Solver::Get_SpeciesDiffusion_TimeStep(Mesher MESH){
int i, j, k, n;

    for (n = 1; n <= 1; n++){

        for(i = Ix[Rango]; i < Fx[Rango]; i++){	
		    for (j = MESH.NY_ColumnMP[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + Halo - Ix[Rango]][1]; j++){
			    for(k = 0; k < NZ; k++){

				    // Mass Diffusion (X Direction)
				    DiffusiveDeltaT += ((CourantFactor * pow(MESH.DeltasMP[LP(i,j,k,0)], 2.0)) / (Species[n].Lewis * Species[n].D_ab[LP(i,j,k,0)] + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * pow(MESH.DeltasMP[LP(i,j,k,0)], 2.0)) / (Species[n].Lewis * Species[n].D_ab[LP(i,j,k,0)] + 1e-10) <= DiffusiveDeltaT);
                  
				    // Mass Diffusion (Y Direction)
				    DiffusiveDeltaT += ((CourantFactor * pow(MESH.DeltasMP[LP(i,j,k,1)], 2.0)) / (Species[n].Lewis * Species[n].D_ab[LP(i,j,k,0)] + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * pow(MESH.DeltasMP[LP(i,j,k,1)], 2.0)) / (Species[n].Lewis * Species[n].D_ab[LP(i,j,k,0)] + 1e-10) <= DiffusiveDeltaT);

				    // Mass Diffusion (Z Direction)
				    DiffusiveDeltaT += ((CourantFactor * pow(MESH.DeltasMP[LP(i,j,k,2)], 2.0)) / (Species[n].Lewis * Species[n].D_ab[LP(i,j,k,0)] + 1e-10) - DiffusiveDeltaT) * ((CourantFactor * pow(MESH.DeltasMP[LP(i,j,k,2)], 2.0)) / (Species[n].Lewis * Species[n].D_ab[LP(i,j,k,0)] + 1e-10) <= DiffusiveDeltaT);

			    }
		    }
	    }

    }
	
	
	MPI_Allreduce(&DiffusiveDeltaT, &DiffusiveDeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

}

// Function to Calculate the Species Equation contributions
void Solver::Get_SpeciesContributions(Mesher MESH){
int i, j, k, n;

    for (n = 0; n < N_Species - 1; n++){

        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
		        for(k = 0; k < NZ; k++){
				    Species[n].ContributionPres[LP(i,j,k,0)] = Species[n].Diffusive[LP(i,j,k,0)] - Species[n].Convective[LP(i,j,k,0)];// + Species[n].wk[LP(i,j,k,0)];
			    }
		    }
	    }

    }

}

// Function to Calculate the New Mass Fraction of the Species
void Solver::Get_Species_MassFraction(Mesher MESH){
int i, j, k, n;

    for (n = 0; n < N_Species - 1; n++){

        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
		        for(k = 0; k < NZ; k++){
				    Species[n].Y_Fut[LP(i,j,k,0)] = Species[n].Y_Pres[LP(i,j,k,0)] + DeltaT * (1.50 * Species[n].ContributionPres[LP(i,j,k,0)] - 0.50 * Species[n].ContributionPast[LP(i,j,k,0)]);
			    }
		    }
	    }

    }
    
}

// Function to set the last species concentration to be (1 - sum(C))
void Solver::Get_Species_MassConservation(Mesher MESH){
int i, j, k, n;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
            for (k = 0; k < NZ; k++){

                Species[N_Species-1].Y_Pres[LP(i,j,k,0)] = 1.0;

                for (n = 0; n < N_Species - 1; n++){
                    Species[N_Species-1].Y_Pres[LP(i,j,k,0)] += - Species[n].Y_Fut[LP(i,j,k,0)];
                }
                
            }
        }
    }

}

// Function to update the mass fraction fields of the species
void Solver::Get_Species_Update(Mesher MESH){
int i, j, k, n;

    for (n = 0; n < N_Species - 1; n++){

        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
                for (k = 0; k < NZ; k++){
                    Species[n].Y_Pres[LP(i,j,k,0)] = Species[n].Y_Fut[LP(i,j,k,0)];
                    Species[n].ContributionPast[LP(i,j,k,0)] = Species[n].ContributionPres[LP(i,j,k,0)];
                }
            }
        }

    }

}