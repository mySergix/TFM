//------------------------------------------------------------------------------------------------//
//                       CPP FILE FOR SPECIES UTILITIES FUNCTIONS                                 //
//------------------------------------------------------------------------------------------------//

// Function to set the initial fields
void Solver::Get_Species_InitialConditions(){
int i, j, k, n;

    // Species n = 0 to n = N - 1
    for (n = 0; n < N_Species - 1; n++){
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
		    for(j = 0; j < NY; j++){
			    for(k = 0; k < NZ; k++){	
				    Species[n].Y_Pres[LP(i,j,k,0)] = Co;
				    Species[n].Y_Fut[LP(i,j,k,0)] = Co;

                    Species[n].ContributionPres[LP(i,j,k,0)] = 0.0;
                    Species[n].ContributionPast[LP(i,j,k,0)] = 0.0;
			    }
		    }
	    }
    }
		
    // Species n = N
    for(i = Ix[Rango]; i < Fx[Rango]; i++){
		for(j = 0; j < NY; j++){
			for(k = 0; k < NZ; k++){	
                Species[N_Species-1].Y_Pres[LP(i,j,k,0)] = 1.0;
                for (n = 0; n < N_Species - 1; n++){
                    Species[N_Species-1].Y_Pres[LP(i,j,k,0)] += - Species[n].Y_Pres[LP(i,j,k,0)];
                }
			}
		}
	}

}

// Function to Calculate the Boussinesq Concetration Buoyancy Term
void Solver::Get_Boussinesq_MassFraction(Mesher MESH, double *V_Matrix){
int i, j, k, n;
double Mass_Fraction;

    // Species 0 -> Water Vapour
    n = 0;
    
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 1; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Mass_Fraction = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                V.Boussinesq_Y[LV(i,j,k,0)] = Beta_Y * V.Gravity * (Mass_Fraction - Co);
            }
        }
    }

}

// Function to Calculate the Species Equation contributions
void Solver::Get_SpeciesContributions(){
int i, j, k, n;

    for (n = 0; n < N_Species - 1; n++){
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = 0; j < NY; j++){
		        for(k = 0; k < NZ; k++){
				    Species[n].ContributionPres[LP(i,j,k,0)] = Species[n].Diffusive[LP(i,j,k,0)] - Species[n].Convective[LP(i,j,k,0)];
			    }
		    }
	    }
    }

}

// Function to Calculate the New Mass Fraction of the Species
void Solver::Get_Species_MassFraction(){
int i, j, k, n;

    for (n = 0; n < N_Species - 1; n++){
        for(i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = 0; j < NY; j++){
		        for(k = 0; k < NZ; k++){
				    Species[n].Y_Fut[LP(i,j,k,0)] = Species[n].Y_Pres[LP(i,j,k,0)] + DeltaT * (1.50 * Species[n].ContributionPres[LP(i,j,k,0)] - 0.50 * Species[n].ContributionPast[LP(i,j,k,0)]);
			    }
		    }
	    }
    }
    
}

// Function to set the last species concentration to be (1 - sum(C))
void Solver::Get_Species_MassConservation(){
int i, j, k, n;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){
                Species[N_Species-1].Y_Pres[LP(i,j,k,0)] = 1.0;
                for (n = 0; n < N_Species - 1; n++){
                    Species[N_Species-1].Y_Pres[LP(i,j,k,0)] += - Species[n].Y_Fut[LP(i,j,k,0)];
                }
                
            }
        }
    }

}

