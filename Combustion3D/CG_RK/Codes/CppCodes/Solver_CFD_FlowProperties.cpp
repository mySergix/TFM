//------------------------------------------------------------------------------------------------//
//                        CPP FILE FOR CFD SOLVER CLASS FLOW PROPERTIES                           //
//------------------------------------------------------------------------------------------------//

// Function to calculate the dynamic viscosity (mu) of each control volume
void Solver::Get_DynamicViscosity(Mesher MESH){
int i, j, k;

	for(i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
			for(k = 0; k < NZ; k++){	
				mu_visc[LP(i,j,k,0)] = JANAF_DynViscosity(T.Pres[LP(i,j,k,0)], i, j, k);
			}
		}
	}

}

// Function to calculate the thermal conduction (Lambda) of each control volume
void Solver::Get_ThermalConductivity(Mesher MESH){
int i, j, k;

	for(i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
			for(k = 0; k < NZ; k++){	
				K_Thermal[LP(i,j,k,0)] = JANAF_ThermalCond(T.Pres[LP(i,j,k,0)], i, j, k);
			}
		}
	}

}

// Function to calculate the Cp Heat Value of ech control volume
void Solver::Get_CpHeat(Mesher MESH){
int i, j, k;

	for(i = Ix[Rango] - 1; i < Fx[Rango] + 1; i++){
		for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
			for(k = 0; k < NZ; k++){	
				Cp_Heat[LP(i,j,k,0)] = JANAF_CpHeat(T.Pres[LP(i,j,k,0)], i, j, k);
			}
		}
	}

}