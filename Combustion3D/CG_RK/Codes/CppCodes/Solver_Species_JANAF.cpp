//------------------------------------------------------------------------------------------------//
//                             CPP FILE FOR JANAF PROPERTIES CALCULATIONS                         //
//------------------------------------------------------------------------------------------------//

// Function to calculate the CP of a certain species
double Solver::JANAF_CpSpecie(double T, int n){
double Cp_Specie;

    if (T > 200.0 && T <= 1000.0){
        Cp_Specie = (R_ideal / Species[n].Wmolar) * (Species[n].Cp_coeff[0] + Species[n].Cp_coeff[1] * T + Species[n].Cp_coeff[2] * pow(T, 2) + Species[n].Cp_coeff[3] * pow(T, 3) + Species[n].Cp_coeff[4] * pow(T, 4));
    }
    else if (T > 1000.0 && T < 5000.0){
        Cp_Specie = (R_ideal / Species[n].Wmolar) * (Species[n].Cp_coeff[5] + Species[n].Cp_coeff[6] * T + Species[n].Cp_coeff[7] * pow(T, 2) + Species[n].Cp_coeff[8] * pow(T, 3) + Species[n].Cp_coeff[9] * pow(T, 4));
    }
    
    return Cp_Specie;

}

// Function to calculate the thermal conductivity of a certain species
double Solver::JANAF_LambdaSpecies(double T, int n){
double lambda;

    if (T > 200.0 && T <= 1000.0){
        lambda = exp(Species[n].lambda_coeff[0] + Species[n].lambda_coeff[1] * log(T) + Species[n].lambda_coeff[2] * pow(log(T), 2) + Species[n].lambda_coeff[3] * pow(log(T), 3));
    }
    else if (T > 1000.0 && T < 5000.0){
        lambda = exp(Species[n].lambda_coeff[4] + Species[n].lambda_coeff[5] * log(T) + Species[n].lambda_coeff[6] * pow(log(T), 2) + Species[n].lambda_coeff[7] * pow(log(T), 3));
    }

    return lambda;

}

// Function to calculate the Heat Capacity coefficient of the mixture in a control volume
double Solver::JANAF_CpHeat(double T, int i, int j, int k){
int n;
double Cp_specie;
double Cp_global = 0.0;

    if (T > 200.0 && T <= 1000.0){
        for (n = 0; n < N_Species; n++){
            Cp_specie = (R_ideal / Species[n].Wmolar) * (Species[n].Cp_coeff[0] + Species[n].Cp_coeff[1] * T + Species[n].Cp_coeff[2] * pow(T, 2) + Species[n].Cp_coeff[3] * pow(T, 3) + Species[n].Cp_coeff[4] * pow(T, 4));
            Cp_global += Species[n].Y_Pres[LP(i,j,k,0)] * Cp_specie;
        }
    }
    else if (T > 1000.0 && T < 5000.0){
        for (n = 0; n < N_Species; n++){
            Cp_specie = (R_ideal / Species[n].Wmolar) * (Species[n].Cp_coeff[5] + Species[n].Cp_coeff[6] * T + Species[n].Cp_coeff[7] * pow(T, 2) + Species[n].Cp_coeff[8] * pow(T, 3) + Species[n].Cp_coeff[9] * pow(T, 4));
            Cp_global += Species[n].Y_Pres[LP(i,j,k,0)] * Cp_specie;
        }
    }
    
    return Cp_global;

}

// Function to calculate the enthalpy of a specie for a certain Temperature
double Solver::JANAF_AbsEnthalpy_Specie(int n, double T){
double H;

    if (T > 200.0 && T <= 1000.0){
        H = (R_ideal / Species[n].Wmolar) * T * (Species[n].h_coeff[0] + (Species[n].h_coeff[1]/2) * T + (Species[n].h_coeff[2]/3) * pow(T, 2) + (Species[n].h_coeff[3]/4) * pow(T, 3) + (Species[n].h_coeff[4]/5) * pow(T, 4) + (Species[n].h_coeff[5]/T));
    }
    else if (T > 1000.0 && T < 5000.0){
        H = (R_ideal / Species[n].Wmolar) * T * (Species[n].h_coeff[6] + (Species[n].h_coeff[7]/2) * T + (Species[n].h_coeff[8]/3) * pow(T, 2) + (Species[n].h_coeff[9]/4) * pow(T, 3) + (Species[n].h_coeff[10]/5) * pow(T, 4) + (Species[n].h_coeff[11]/T));
    }
    
    return H;
    
}

// Function to calculate the absolute enthalpy of the mixture in a control volume
double Solver::JANAF_AbsEnthalpy_Specie_Mix(double T, int i, int j, int k){
int n;
double AbsEnth;
double h = 0.0;

    if (T > 200.0 && T <= 1000.0){
        for (n = 0; n < N_Species; n++){
            AbsEnth = (R_ideal / Species[n].Wmolar) * T * (Species[n].h_coeff[0] + (Species[n].h_coeff[1]/2) * T + (Species[n].h_coeff[2]/3) * pow(T, 2) + (Species[n].h_coeff[3]/4) * pow(T, 3) + (Species[n].h_coeff[4]/5) * pow(T, 4) + (Species[n].h_coeff[5]/T));
            h += Species[n].Y_Pres[LP(i,j,k,0)] * AbsEnth;
        }
    }
    else if (T > 1000.0 && T < 5000.0){
        for (n = 0; n < N_Species; n++){
            AbsEnth = (R_ideal / Species[n].Wmolar) * T * (Species[n].h_coeff[6] + (Species[n].h_coeff[7]/2) * T + (Species[n].h_coeff[8]/3) * pow(T, 2) + (Species[n].h_coeff[9]/4) * pow(T, 3) + (Species[n].h_coeff[10]/5) * pow(T, 4) + (Species[n].h_coeff[11]/T));
            h += Species[n].Y_Pres[LP(i,j,k,0)] * AbsEnth;
        }
    }
    
    return h;

}

// Function to calculate the dynamic viscosity of the mixture in a control volume
double Solver::JANAF_DynViscosity(double T, int i, int j, int k){
int n;
double DynVisc;
double mu = 0.0;

    if (T > 200.0 && T <= 1000.0){
        for (n = 0; n < N_Species; n++){
            DynVisc = exp(Species[n].mu_coeff[0] + Species[n].mu_coeff[1] * log(T) + Species[n].mu_coeff[2] * pow(log(T), 2) + Species[n].mu_coeff[3] * pow(log(T), 3));
            mu += Species[n].Y_Pres[LP(i,j,k,0)] * DynVisc;
        }
    }
    else if (T > 1000.0 && T < 5000.0){
        for (n = 0; n < N_Species; n++){
            DynVisc = exp(Species[n].mu_coeff[4] + Species[n].mu_coeff[5] * log(T) + Species[n].mu_coeff[6] * pow(log(T), 2) + Species[n].mu_coeff[7] * pow(log(T), 3));
            mu += Species[n].Y_Pres[LP(i,j,k,0)] * DynVisc;
        }
    }

    return mu;

}

// Function to calculate the thermal conductivity of the mixture in a control volume
double Solver::JANAF_ThermalCond(double T, int i, int j, int k){
int n;
double ThermalCond;
double lambda = 0.0;

    if (T > 200.0 && T <= 1000.0){
        for (n = 0; n < N_Species; n++){
            ThermalCond = exp(Species[n].lambda_coeff[0] + Species[n].lambda_coeff[1] * log(T) + Species[n].lambda_coeff[2] * pow(log(T), 2) + Species[n].lambda_coeff[3] * pow(log(T), 3));
            lambda += Species[n].Y_Pres[LP(i,j,k,0)] * ThermalCond;
        }
    }
    else if (T > 1000.0 && T < 5000.0){
        for (n = 0; n < N_Species; n++){
            ThermalCond = exp(Species[n].lambda_coeff[4] + Species[n].lambda_coeff[5] * log(T) + Species[n].lambda_coeff[6] * pow(log(T), 2) + Species[n].lambda_coeff[7] * pow(log(T), 3));
            lambda += Species[n].Y_Pres[LP(i,j,k,0)] * ThermalCond;
        }
    }

    return lambda;

}