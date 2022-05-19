//------------------------------------------------------------------------------------------------//
//                       CPP FILE FOR SPECIES CHEMISTRY FUNCTIONS                                 //
//------------------------------------------------------------------------------------------------//

// Function to get all the reaction kinetics data
void Solver::Get_ReactionsData_OneStepCH4(){

    Reactions.A[0] = 4.2e15;
    Reactions.Beta[0] = 0.0;
    Reactions.EA[0] = 16900.0;

    Reactions.DeltaS[0] = 0.0;
    Reactions.DeltaH[0] = 0.0;

}

// Function to calculate the reaction rate of progress of each species (4 Step - CH4)
void Solver::Get_Reactions_OneStepCH4(Mesher MESH){
int i, j, k, n, r;
double Q1;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0]; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1]; j++){
            for (k = 0; k < NZ; k++){
                Reactions.Kf[0] = Reactions.A[0] * exp(- Reactions.EA[0] / 298.15);
                Q1 = Reactions.Kf[0] * pow(Rho * Species[0].Y_Pres[LP(i,j,k,0)] / Species[0].Wmolar, 2.8) * pow(Rho * Species[1].Y_Pres[LP(i,j,k,0)] / Species[1].Wmolar, 1.2);
                
                // Species 0 - CH4
                Species[0].wk[LP(i,j,k,0)] = -1.0 * Q1 * Species[0].Wmolar;

                // Species 1 - O2
                Species[1].wk[LP(i,j,k,0)] = -2.0 * Q1 * Species[1].Wmolar;

                // Species 2 - CO2
                Species[2].wk[LP(i,j,k,0)] = 1.0 * Q1 * Species[2].Wmolar;

                // Species 3 - H2O
                Species[3].wk[LP(i,j,k,0)] = 2.0 * Q1 * Species[3].Wmolar;

            }
        }
    }

}

// Function to calculate the reaction rate of progress of each species (4 Step - CH4)
void Solver::Get_Reactions_FourStepCH4(Mesher MESH){
int i, j, k, n, r;
double Q1, Q2, Q3, Q4;

    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (j = 0; j < NY; j++){
            for (k = 0; k < NZ; k++){

                // Reaction rate coefficients calculation
                for (r = 0; r < N_Reactions; r++){
                    Reactions.Kf[r] = Reactions.A[r] * pow(T.Pres[LP(i,j,k,0)], Reactions.Beta[r]) * exp(- Reactions.EA[r] / (R_ideal * T.Pres[LP(i,j,k,0)]));
                }

                // Reaction 3 (Reversible)
                Reactions.Kr[2] = Reactions.Kf[2] / ((P.Pres[LP(i,j,k,0)] / (R_ideal * T.Pres[LP(i,j,k,0)])) * exp((Reactions.DeltaS[2] / R_ideal) - Reactions.DeltaH[2] / (R_ideal * T.Pres[LP(i,j,k,0)])));

                // Reaction 4 (Reversible)
                Reactions.Kr[3] = Reactions.Kf[3] / ((P.Pres[LP(i,j,k,0)] / (R_ideal * T.Pres[LP(i,j,k,0)])) * exp((Reactions.DeltaS[3] / R_ideal) - Reactions.DeltaH[3] / (R_ideal * T.Pres[LP(i,j,k,0)])));

                
                // Reactions rates calculation
                Q1 = Reactions.Kf[0] * pow(Rho * Species[0].Y_Pres[LP(i,j,k,0)] / Species[0].Wmolar, 0.50) * pow(Rho * Species[1].Y_Pres[LP(i,j,k,0)] / Species[1].Wmolar, 1.25);
                
                Q2 = Reactions.Kf[1] * pow(Rho * Species[0].Y_Pres[LP(i,j,k,0)] / Species[0].Wmolar, 1.0) * pow(Rho * Species[4].Y_Pres[LP(i,j,k,0)] / Species[4].Wmolar, 1.0);
                
                Q3 = Reactions.Kf[2] * pow(Rho * Species[1].Y_Pres[LP(i,j,k,0)] / Species[1].Wmolar, 2.25) * pow(Rho * Species[5].Y_Pres[LP(i,j,k,0)] / Species[5].Wmolar, 0.50) * pow(Rho * Species[4].Y_Pres[LP(i,j,k,0)] / Species[4].Wmolar, -1.0)
                   - Reactions.Kr[2] * pow(Rho * Species[1].Y_Pres[LP(i,j,k,0)] / Species[1].Wmolar, 1.75) * pow(Rho * Species[5].Y_Pres[LP(i,j,k,0)] / Species[5].Wmolar, -0.50);

                Q4 = Reactions.Kf[3] * pow(Rho * Species[2].Y_Pres[LP(i,j,k,0)] / Species[2].Wmolar, 1.0) * pow(Rho * Species[4].Y_Pres[LP(i,j,k,0)] / Species[4].Wmolar, 1.0)
                   - Reactions.Kr[3] * pow(Rho * Species[3].Y_Pres[LP(i,j,k,0)] / Species[3].Wmolar, 1.0) * pow(Rho * Species[4].Y_Pres[LP(i,j,k,0)] / Species[4].Wmolar, 1.0);
            
                // Species Production/Destruction Calculation
                
                    // CH4
                    Species[0].wk[LP(i,j,k,0)] = -1.0 * Q1 - 1.0 * Q2;

                    // O2
                    Species[1].wk[LP(i,j,k,0)] = -0.50 * Q1 - 0.50 * Q3;

                    // CO
                    Species[2].wk[LP(i,j,k,0)] = 1.0 * Q1 + 1.0 * Q2 - 1.0 * Q4;

                    // CO2
                    Species[3].wk[LP(i,j,k,0)] = 1.0 * Q4;

                    // H2O
                    Species[4].wk[LP(i,j,k,0)] = -1.0 * Q2 + 1.0 * Q3 - 1.0 * Q4;

                    // H2
                    Species[5].wk[LP(i,j,k,0)] = 2.0 * Q1 + 3.0 * Q2 - 1.0 * Q3 + 1.0 * Q4;
            
            }
        }
    }

}