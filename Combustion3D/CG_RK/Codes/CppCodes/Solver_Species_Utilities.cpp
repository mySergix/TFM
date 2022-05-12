//------------------------------------------------------------------------------------------------//
//                       CPP FILE FOR SPECIES UTILITIES FUNCTIONS                                 //
//------------------------------------------------------------------------------------------------//

// Function to set the initial fields
void Solver::Get_Species_InitialConditions(Mesher MESH){
int i, j, k, n;

    for (n = 0; n < N_Species; n++){

        // Fluid Domain
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
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
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
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
		
}

// Function to Calculate the Species Equation contributions
void Solver::Get_SpeciesContributions(Mesher MESH){
int i, j, k, n;

    for (n = 0; n < N_Species - 1; n++){

        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
		        for(k = 0; k < NZ; k++){
				    Species[n].ContributionPres[LP(i,j,k,0)] = Species[n].Diffusive[LP(i,j,k,0)] - Species[n].Convective[LP(i,j,k,0)];
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

