//------------------------------------------------------------------------------------------------//
//                        CPP FILE FOR NODAL COORDINATES DISCRETIZATION                           //
//------------------------------------------------------------------------------------------------//

// Nodal coordinates discretization
void Mesher::Get_LocalMeshes(){
int i, j, k;
double I, J, K;
double nx;
double ny;
double nz;

	// Coordinates X

	//Staggered U Mesh	
	for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo + 1; i++){
		I = i;
		for (j = - Halo; j < NY + Halo; j++){
			J = j;
			for (k = - Halo; k < NZ + Halo; k++){

				// Section 1 (X Direction)
				if (i <= NX_1){
					nx = NX_1;
					// Regular
					if (OptionX_1 == 1){
						MU[LU(i,j,k,0)] = (I - Ix[Rango]) * (X_1 / nx);
					}
					// Right Sided - Hyperbolic Tangent
					else if (OptionX_1 == 2){
						MU[LU(i,j,k,0)] = (X_1 / 2.0) * (tanh(SFX_1 * ((I - Ix[Rango]) / nx)) / tanh(SFX_1));
					}
					
				}
				// Section 2 (X Direction)
				else if (i > NX_1 && i <= NX_1 + NX_2){
					nx = NX_2;
					// Regular
					if (OptionX_2 == 1){
						MU[LU(i,j,k,0)] = X_1 + (I - Ix[Rango]) * (X_2 / nx);
					}
					// Right Sided - Hyperbolic Tangent
					else if (OptionX_2 == 2){
						MU[LU(i,j,k,0)] = X_1 + (X_2 / 2.0) * (tanh(SFX_2 * ((I - Ix[Rango]) / nx)) / tanh(SFX_2));
					}

				}
				// Section 3 (X Direction)
				else if (i > NX_1 + NX_2){
					nx = NX_3;
					// Regular
					if (OptionX_3 == 1){
						MU[LU(i,j,k,0)] = X_1 + X_2 + (I - Ix[Rango]) * (X_3 / nx);
					}
					// Right Sided - Hyperbolic Tangent
					else if (OptionX_3 == 2){
						MU[LU(i,j,k,0)] = X_1 + X_2 + (X_3 / 2.0) * (tanh(SFX_3 * ((I - Ix[Rango]) / nx)) / tanh(SFX_3));
					}
					// Left Sided - Hyperbolic Tangent
					else if (OptionX_3 == 3){
						MU[LU(i,j,k,0)] = X_1 + X_2 + (X_3 / 2.0) * (tanh(SFX_3 * ((nx - (I - Ix[Rango])) / nx)) / tanh(SFX_3));
					}
					// Centered - Hyperbolic Tangent
					else if (OptionX_3 == 4){
						MU[LU(i,j,k,0)] = X_1 + X_2 + (X_3 / 2.0) * (tanh(SFX_3 * (((2.0 * (I - Ix[Rango])) - nx) / nx)) / tanh(SFX_3));
					}
				}

			}
		}
	}

	// Collocated Pressure mesh
	for(i = Ix[Rango] - HP; i < Fx[Rango] + HP; i++){
		for(j = - HP; j < NY + HP; j++){
			for(k = - HP; k < NZ + HP; k++){
				MP[LP(i,j,k,0)] = 0.50*(MU[LU(i,j,k,0)] + MU[LU(i+1,j,k,0)]);
			}
		}
	}

	// Staggered V mesh
	for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
		for(k = - Halo; k < NZ + Halo; k++){
			MV[LV(i,NY + Halo,k,0)] = 0.50*(MU[LU(i,NY - 1 + Halo,k,0)] + MU[LU(i + 1,NY-1 + Halo,k,0)]);
			for(j = - Halo; j < NY + Halo; j++){
				MV[LV(i,j,k,0)] = 0.50*(MU[LU(i,j,k,0)] + MU[LU(i+1,j,k,0)]);
			}
		}
	}

	// Staggered W mesh
	for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
		for(j = - Halo; j < NY + Halo; j++){
			MW[LW(i,j,NZ + Halo,0)] = 0.50*(MU[LU(i,j,NZ-1 + Halo,0)] + MU[LU(i+1,j,NZ-1 + Halo,0)]);
			for(k = - Halo; k < NZ + Halo; k++){
				MW[LW(i,j,k,0)] = 0.50*(MU[LU(i,j,k,0)] + MU[LU(i+1,j,k,0)]);
			}
		}
	}

	
	// Coordinates Y

	// Staggered V Mesh
	if (Problema == "Premixed"){

		for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
			I = i;
			for (j = - Halo; j < NY + Halo + 1; j++){
				J = j;
				for (k = - Halo; k < NZ + Halo; k++){

					// Section 1 (Y Direction)
					if (j <= NY_1){
						ny = NY_1;
						// Regular
						if (OptionY_1 == 1){
							MV[LV(i,j,k,1)] = J * (Y_1 / ny);
						}
						// Top Sided - Hyperbolic Tangent
						else if (OptionY_1 == 2){
							MV[LV(i,j,k,1)] = (Y_1 / 2.0) * (tanh(SFY_1 * (J / ny)) / tanh(SFY_1));
						}
					
					}
					// Section 2 (Y Direction)
					else if (j > NY_1 && j <= NY_1 + NY_2){
						ny = NY_2;
						// Regular
						if (OptionY_2 == 1){
							MV[LV(i,j,k,1)] = Y_1 + J * (Y_2 / ny);
						}
						// Top Sided - Hyperbolic Tangent
						else if (OptionY_2 == 2){
							MV[LV(i,j,k,1)] = Y_1 + (Y_2 / 2.0) * (tanh(SFY_2 * (J / ny)) / tanh(SFY_2));
						}

					}
					// Section 3 (Y Direction)
					else if (j > NY_1 + NY_2){
						ny = NY_3;
						// Regular
						if (OptionY_3 == 1){
							MV[LV(i,j,k,1)] = Y_1 + Y_2 + J * (Y_3 / ny);
						}
						// Right Sided - Hyperbolic Tangent
						else if (OptionY_3 == 2){
							MV[LV(i,j,k,1)] = Y_1 + Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * (J / ny)) / tanh(SFY_3));
						}
						// Left Sided - Hyperbolic Tangent
						else if (OptionY_3 == 3){
							MV[LV(i,j,k,1)] = Y_1 + Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * ((ny - J) / ny)) / tanh(SFY_3));
						}
						// Centered - Hyperbolic Tangent
						else if (OptionY_3 == 4){
							MV[LV(i,j,k,1)] = Y_1 + Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * (((2.0 * J) - ny) / ny)) / tanh(SFY_3));
						}
					}

				}
			}
		}

	}
	else if (Problema == "NonPremixed"){

		for (i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
			I = i;
			for (j = - Halo; j < NY + Halo + 1; j++){
				J = j;
				for (k = - Halo; k < NZ + Halo; k++){

					// Section 2 (Y Direction)
					if (j <= NY_2){
						ny = NY_2;
						// Regular
						if (OptionY_2 == 1){
							MV[LV(i,j,k,1)] = J * (Y_2 / ny);
						}
						// Top Sided - Hyperbolic Tangent
						else if (OptionY_2 == 2){
							MV[LV(i,j,k,1)] = (Y_2 / 2.0) * (tanh(SFY_2 * (J / ny)) / tanh(SFY_2));
						}

					}
					// Section 3 (Y Direction)
					else if (j > NY_2){
						ny = NY_3;
						// Regular
						if (OptionY_3 == 1){
							MV[LV(i,j,k,1)] = Y_2 + J * (Y_3 / ny);
						}
						// Right Sided - Hyperbolic Tangent
						else if (OptionY_3 == 2){
							MV[LV(i,j,k,1)] = Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * (J / ny)) / tanh(SFY_3));
						}
						// Left Sided - Hyperbolic Tangent
						else if (OptionY_3 == 3){
							MV[LV(i,j,k,1)] = Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * ((ny - J) / ny)) / tanh(SFY_3));
						}
						// Centered - Hyperbolic Tangent
						else if (OptionY_3 == 4){
							MV[LV(i,j,k,1)] = Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * (((2.0 * J) - ny) / ny)) / tanh(SFY_3));
						}
					}

				}
			}
		}

	}

	// Collocated Pressure mesh
	for(i = Ix[Rango] - HP; i < Fx[Rango] + HP; i++){
		for(j = - HP; j < NY + HP; j++){
			for(k = - HP; k < NZ + HP; k++){
				MP[LP(i,j,k,1)] = 0.50*(MV[LV(i,j,k,1)] + MV[LV(i,j+1,k,1)]);
			}
		}
	}

	// Staggered U mesh
	for(k = - Halo; k < NZ + Halo; k++){
		for(j = - Halo; j < NY + Halo; j++){
			MU[LU(Fx[Rango] + Halo,j,k,1)] = 0.50*(MV[LV(Fx[Rango]-1 + Halo,j,k,1)] + MV[LV(Fx[Rango]-1 + Halo,j+1,k,1)]);
			for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
				MU[LU(i,j,k,1)] = 0.50*(MV[LV(i,j,k,1)] + MV[LV(i,j+1,k,1)]);
			}
		}
	}

	// Staggered W mesh
	for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
		for(j = - Halo; j < NY + Halo; j++){
			MW[LW(i,j,NZ + Halo,1)] = 0.50*(MV[LV(i,j,NZ-1 + Halo,1)] + MV[LV(i,j+1,NZ-1 + Halo,1)]);
			for(k = - Halo; k < NZ + Halo; k++){
				MW[LW(i,j,k,1)] = 0.50*(MV[LV(i,j,k,1)] + MV[LV(i,j+1,k,1)]);
			}
		}
	}   


	// Coordinates Z

	// Staggered W mesh (Regular)
	nz = NZ;
	for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
		for(j = - Halo; j < NY + Halo; j++){
			for(k = - Halo; k < NZ + Halo + 1; k++){
				K = k;
				MW[LW(i,j,k,2)] = K * (Zdomain / nz);
			}
		}
	}

	// Collocated Pressure mesh
	for(i = Ix[Rango] - HP; i < Fx[Rango] + HP; i++){
		for(j = - Halo; j < NY + Halo; j++){
			for(k = - Halo; k < NZ + Halo; k++){
				MP[LP(i,j,k,2)] = 0.50*(MW[LW(i,j,k,2)] + MW[LW(i,j,k+1,2)]);
			}
		}
	}

	// Staggered U mesh
	for(k = - Halo; k < NZ + Halo; k++){
		for(j = - Halo; j < NY + Halo; j++){
			MU[LU(Fx[Rango] + Halo,j,k,2)] = 0.50*(MW[LW(Fx[Rango] - 1 + Halo,j,k,2)] + MW[LW(Fx[Rango] - 1 + Halo,j,k+1,2)]);
			for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
				MU[LU(i,j,k,2)] = 0.50*(MW[LW(i,j,k,2)] + MW[LW(i,j,k+1,2)]);
			}
		}
	}

	// Staggered V mesh
	for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
		for(k = - Halo; k < NZ + Halo; k++){
			MV[LV(i,NY + Halo,k,2)] = 0.50*(MW[LW(i,NY-1 + Halo,k,2)] + MW[LW(i,NY-1 + Halo,k+1,2)]);
			for(j = - Halo; j < NY + Halo; j++){
				MV[LV(i,j,k,2)] = 0.50*(MW[LW(i,j,k,2)] + MW[LW(i,j,k+1,2)]);
			}
		}
	}

}

// Function to calculate the coordinates of the global mesh
void Mesher::Get_GlobalMesh(){
int i, j, k;
double I, J, K;
double nx = NX;
double ny = NY;
double nz = NZ;

	// Coordinates X

	//Staggered U Mesh	
	for (i = 0; i < NX + 1; i++){
		I = i;
		for (j = 0; j < NY; j++){
			J = j;
			for (k = 0; k < NZ; k++){

				// Section 1 (X Direction)
				if (i <= NX_1){
					nx = NX_1;
					// Regular
					if (OptionX_1 == 1){
						GlobalMeshU[GU(i,j,k,0)] = I * (X_1 / nx);
					}
					// Right Sided - Hyperbolic Tangent
					else if (OptionX_1 == 2){
						GlobalMeshU[GU(i,j,k,0)] = (X_1 / 2.0) * (tanh(SFX_1 * (I / nx)) / tanh(SFX_1));
					}
					
				}
				// Section 2 (X Direction)
				else if (i >= NX_1 && i <= NX_1 + NX_2){
					nx = NX_2;
					// Regular
					if (OptionX_2 == 1){
						GlobalMeshU[GU(i,j,k,0)] = X_1 + (I - NX_1) * (X_2 / nx);
					}
					// Right Sided - Hyperbolic Tangent
					else if (OptionX_2 == 2){
						GlobalMeshU[GU(i,j,k,0)] = X_1 + (X_2 / 2.0) * (tanh(SFX_2 * (I / nx)) / tanh(SFX_2));
					}

				}
				// Section 3 (X Direction)
				else if (i >= NX_1 + NX_2){
					nx = NX_3;
					// Regular
					if (OptionX_3 == 1){
						GlobalMeshU[GU(i,j,k,0)] = X_1 + X_2 + (I - (NX_1 + NX_2)) * (X_3 / nx);
					}
					// Right Sided - Hyperbolic Tangent
					else if (OptionX_3 == 2){
						GlobalMeshU[GU(i,j,k,0)] = X_1 + X_2 + (X_3 / 2.0) * (tanh(SFX_3 * (I / nx)) / tanh(SFX_3));
					}
					// Left Sided - Hyperbolic Tangent
					else if (OptionX_3 == 3){
						GlobalMeshU[GU(i,j,k,0)] = X_1 + X_2 + (X_3 / 2.0) * (tanh(SFX_3 * ((nx - I) / nx)) / tanh(SFX_3));
					}
					// Centered - Hyperbolic Tangent
					else if (OptionX_3 == 4){
						GlobalMeshU[GU(i,j,k,0)] = X_1 + X_2 + (X_3 / 2.0) * (tanh(SFX_3 * (((2.0 * I) - nx) / nx)) / tanh(SFX_3));
					}
				}

			}
		}
	}

	// Collocated Pressure mesh
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){
			for (k = 0; k < NZ; k++){
				GlobalMeshP[GP(i,j,k,0)] = 0.50 * (GlobalMeshU[GU(i,j,k,0)] + GlobalMeshU[GU(i+1,j,k,0)]);
			}
		}
	}

	// Staggered V mesh
	for (i = 0; i < NX; i++){
		for (k = 0; k < NZ; k++){

			// Top Part
			GlobalMeshV[GV(i,NY,k,0)] = 0.50 * (GlobalMeshU[GU(i,NY-1,k,0)] + GlobalMeshU[GU(i+1,NY-1,k,0)]);

			for (j = 0; j < NY; j++){
				GlobalMeshV[GV(i,j,k,0)] = 0.50 * (GlobalMeshU[GU(i,j,k,0)] + GlobalMeshU[GU(i+1,j,k,0)]);
			}

		}
	}

	// Staggered W mesh
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){

			// There Part
			GlobalMeshW[GW(i,j,NZ,0)] = 0.50 * (GlobalMeshU[GU(i,j,NZ-1,0)] + GlobalMeshU[GU(i+1,j,NZ-1,0)]);

			for (k = 0; k < NZ; k++){
				GlobalMeshW[GW(i,j,k,0)] = 0.50 * (GlobalMeshU[GU(i,j,k,0)] + GlobalMeshU[GU(i+1,j,k,0)]);
			}

		}
	}


	// Coordinates Y

	// Staggered V Mesh
	if (Problema == "Premixed"){

		for (i = 0; i < NX; i++){
			I = i;
			for (j = 0; j < NY + 1; j++){
				J = j;
				for (k = 0; k < NZ; k++){

					// Section 1 (Y Direction)
					if (j <= NY_1){
						ny = NY_1;
						// Regular
						if (OptionY_1 == 1){
							GlobalMeshV[GV(i,j,k,1)] = J * (Y_1 / ny);
						}
						// Top Sided - Hyperbolic Tangent
						else if (OptionY_1 == 2){
							GlobalMeshV[GV(i,j,k,1)] = (Y_1 / 2.0) * (tanh(SFY_1 * (J / ny)) / tanh(SFY_1));
						}
					
					}
					// Section 2 (Y Direction)
					else if (j > NY_1 && j <= NY_1 + NY_2){
						ny = NY_2;
						// Regular
						if (OptionY_2 == 1){
							GlobalMeshV[GV(i,j,k,1)] = Y_1 + (J - NY_1) * (Y_2 / ny);
						}
						// Top Sided - Hyperbolic Tangent
						else if (OptionY_2 == 2){
							GlobalMeshV[GV(i,j,k,1)] = Y_1 + (Y_2 / 2.0) * (tanh(SFY_2 * (J / ny)) / tanh(SFY_2));
						}

					}
					// Section 3 (Y Direction)
					else if (j > NY_1 + NY_2){
						ny = NY_3;
						// Regular
						if (OptionY_3 == 1){
							GlobalMeshV[GV(i,j,k,1)] = Y_1 + Y_2 + (J - (NY_1 + NY_2)) * (Y_3 / ny);
						}
						// Right Sided - Hyperbolic Tangent
						else if (OptionY_3 == 2){
							GlobalMeshV[GV(i,j,k,1)] = Y_1 + Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * (J / ny)) / tanh(SFY_3));
						}
						// Left Sided - Hyperbolic Tangent
						else if (OptionY_3 == 3){
							GlobalMeshV[GV(i,j,k,1)] = Y_1 + Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * ((ny - J) / ny)) / tanh(SFY_3));
						}
						// Centered - Hyperbolic Tangent
						else if (OptionY_3 == 4){
							GlobalMeshV[GV(i,j,k,1)] = Y_1 + Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * (((2.0 * J) - ny) / ny)) / tanh(SFY_3));
						}
					}

				}
			}
		}

	}
	else if (Problema == "NonPremixed"){

		for (i = 0; i < NX; i++){
			I = i;
			for (j = 0; j < NY + 1; j++){
				J = j;
				for (k = 0; k < NZ; k++){

					// Section 2 (Y Direction)
					if (j <= NY_2){
						ny = NY_2;
						// Regular
						if (OptionY_2 == 1){
							MV[LV(i,j,k,1)] = J * (Y_2 / ny);
						}
						// Top Sided - Hyperbolic Tangent
						else if (OptionY_2 == 2){
							MV[LV(i,j,k,1)] = (Y_2 / 2.0) * (tanh(SFY_2 * (J / ny)) / tanh(SFY_2));
						}

					}
					// Section 3 (Y Direction)
					else if (j > NY_2){
						ny = NY_3;
						// Regular
						if (OptionY_3 == 1){
							MV[LV(i,j,k,1)] = Y_2 + J * (Y_3 / ny);
						}
						// Right Sided - Hyperbolic Tangent
						else if (OptionY_3 == 2){
							MV[LV(i,j,k,1)] = Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * (J / ny)) / tanh(SFY_3));
						}
						// Left Sided - Hyperbolic Tangent
						else if (OptionY_3 == 3){
							MV[LV(i,j,k,1)] = Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * ((ny - J) / ny)) / tanh(SFY_3));
						}
						// Centered - Hyperbolic Tangent
						else if (OptionY_3 == 4){
							MV[LV(i,j,k,1)] = Y_2 + (Y_3 / 2.0) * (tanh(SFY_3 * (((2.0 * J) - ny) / ny)) / tanh(SFY_3));
						}
					}

				}
			}
		}

	}

	// Collocated Pressure mesh
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){
			for (k = 0; k < NZ; k++){
				GlobalMeshP[GP(i,j,k,1)] = 0.50 * (GlobalMeshV[GV(i,j,k,1)] + GlobalMeshV[GV(i,j+1,k,1)]);
			}
		}
	}

	// Staggered U mesh
	for (k = 0; k < NZ; k++){
		for (j = 0; j < NY; j++){

			// Right Part
			GlobalMeshU[GU(NX,j,k,1)] = 0.50 * (GlobalMeshV[GV(NX-1,j,k,1)] + GlobalMeshV[GV(NX-1,j+1,k,1)]);

			for (i = 0; i < NX; i++){
				GlobalMeshU[GU(i,j,k,1)] = 0.50 * (GlobalMeshV[GV(i,j,k,1)] + GlobalMeshV[GV(i,j+1,k,1)]);
			}

		}
	}

	// Staggered W mesh
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){

			// There Part
			GlobalMeshW[GW(i,j,NZ,1)] = 0.50 * (GlobalMeshV[GV(i,j,NZ-1,1)] + GlobalMeshV[GV(i,j+1,NZ-1,1)]);

			for (k = 0; k < NZ; k++){
				GlobalMeshW[GW(i,j,k,1)] = 0.50 * (GlobalMeshV[GV(i,j,k,1)] + GlobalMeshV[GV(i,j+1,k,1)]);
			}

		}
	}  
 

	// Coordinates Z

	// Staggered W mesh (Regular)
	nz = NZ;
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){
			for (k = 0; k < NZ + 1; k++){
				K = k;
				GlobalMeshW[GW(i,j,k,2)] = K * (Zdomain / nz);
			}
		}
	}

	// Collocated Pressure mesh
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){
			for (k = 0; k < NZ; k++){
				GlobalMeshP[GP(i,j,k,2)] = 0.50 * (GlobalMeshW[GW(i,j,k,2)] + GlobalMeshW[GW(i,j,k+1,2)]);
			}
		}
	}

	// Staggered U mesh
	for (k = 0; k < NZ; k++){
		for (j = 0; j < NY; j++){

			// Right Part
			GlobalMeshU[GU(NX,j,k,2)] = 0.50 * (GlobalMeshW[GW(NX-1,j,k,2)] + GlobalMeshW[GW(NX-1,j,k+1,2)]);

			for (i = 0; i < NX; i++){
				GlobalMeshU[GU(i,j,k,2)] = 0.50 * (GlobalMeshW[GW(i,j,k,2)] + GlobalMeshW[GW(i,j,k+1,2)]);
			}

		}
	}

	// Staggered V mesh
	for (i = 0; i < NX; i++){
		for (k = 0; k < NZ; k++){

			// Top Part
			GlobalMeshV[GV(i,NY,k,2)] = 0.50 * (GlobalMeshW[GW(i,NY-1,k,2)] + GlobalMeshW[GW(i,NY-1,k+1,2)]);

			for (j = 0; j < NY; j++){
				GlobalMeshV[GV(i,j,k,2)] = 0.50 * (GlobalMeshW[GW(i,j,k,2)] + GlobalMeshW[GW(i,j,k+1,2)]);
			}

		}
	}

}