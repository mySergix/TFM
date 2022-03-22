//------------------------------------------------------------------------------------------------//
//                        CPP FILE FOR NODAL COORDINATES DISCRETIZATION                           //
//------------------------------------------------------------------------------------------------//

// Nodal coordinates discretization
void Mesher::Get_LocalMeshes(){
int i, j, k;
double I, J, K;
double nx = NX;
double ny = NY;
double nz = NZ;

	// Coordinates X

		if(OptionX == 1){ // Regular

			// Staggered U mesh
			for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo + 1; i++){
				for(j = - Halo; j < NY + Halo; j++){
					for(k = - Halo; k < NZ + Halo; k++){
						MU[LU(i,j,k,0)] = i*(Xdominio/nx);
					}
				}
			}

		}
		else if(OptionX == 2){ // Hyperbolic Tangent

			// Staggered U mesh
			for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo + 1; i++){
				for(j = - Halo; j < NY + Halo; j++){
					for(k = - Halo; k < NZ + Halo; k++){
						I = i;
						MU[LU(i,j,k,0)] = (Xdominio/2.0)*(1.0 + (tanh(SFX*((2.0*I - nx)/nx)) + tanh(SFX))/tanh(SFX) - 1.0);
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

		if(OptionY == 1){ // Regular

			// Staggered V mesh
			for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
				for(j = - Halo; j < NY + Halo + 1; j++){
					for(k = - Halo; k < NZ + Halo; k++){
						MV[LV(i,j,k,1)] = j*(Ydominio/ny);
					}
				}
			}

		}
		else if(OptionY == 2){ // Hyperbolic Tangent

			// Staggered V mesh
			for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
				for(j = - Halo; j < NY + Halo + 1; j++){
					for(k = - Halo; k < NZ + Halo; k++){
						J = j;
						MV[LV(i,j,k,1)] = (Ydominio/2.0)*(1.0 + (tanh(SFY*((2.0*J - ny)/ny)) + tanh(SFY))/tanh(SFY) - 1.0);
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

		if(OptionZ == 1){ // Regular

			// Staggered W mesh
			for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
				for(j = - Halo; j < NY + Halo; j++){
					for(k = - Halo; k < NZ + Halo + 1; k++){
						MW[LW(i,j,k,2)] = k*(Zdominio/nz);
					}
				}
			}

		}
		else if(OptionZ == 2){ // Hyperbolic Tangent

			// Staggered W mesh
			for(i = Ix[Rango] - Halo; i < Fx[Rango] + Halo; i++){
				for(j = - Halo; j < NY + Halo; j++){
					for(k = - Halo; k < NZ + Halo + 1; k++){
						K = k;
						MW[LW(i,j,k,2)] = (Zdominio/2.0)*(1.0 + (tanh(SFZ*((2.0*K - nz)/nz)) + tanh(SFZ))/tanh(SFZ) - 1.0);
					}
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

// Function to calculate the coordinates of teh global mesh
void Mesher::Get_GlobalMesh(){
int i, j, k;
double I, J, K;
double nx = NX;
double ny = NY;
double nz = NZ;

	// Coordinates X

	if(OptionX == 1){ // Regular

		// Staggered U mesh
		for (i = 0; i < NX + 1; i++){
			for (j = 0; j < NY; j++){
				for (k = 0; k < NZ; k++){
					GlobalMeshU[GU(i,j,k,0)] = i*(Xdominio/nx);
				}
			}
		}

	}
	else if(OptionX == 2){ // Hyperbolic Tangent

		// Staggered U mesh
		for (i = 0; i < NX + 1; i++){
			for (j = 0; j < NY; j++){
				for (k = 0; k < NZ; k++){
					I = i;
					GlobalMeshU[GU(i,j,k,0)] = (Xdominio/2.0)*(1.0 + (tanh(SFX*((2.0*I - nx)/nx)) + tanh(SFX))/tanh(SFX) - 1.0);
				}
			}
		}

	}

	// Collocated Pressure mesh
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){
			for (k = 0; k < NZ; k++){
				GlobalMeshP[GP(i,j,k,0)] = 0.50*(GlobalMeshU[GU(i,j,k,0)] + GlobalMeshU[GU(i+1,j,k,0)]);
			}
		}
	}

	// Staggered V mesh
	for (i = 0; i < NX; i++){
		for (k = 0; k < NZ; k++){
			GlobalMeshV[GV(i,NY,k,0)] = 0.50*(GlobalMeshU[GU(i,NY-1,k,0)] + GlobalMeshU[GU(i+1,NY-1,k,0)]);
			for (j = 0; j < NY; j++){
				GlobalMeshV[GV(i,j,k,0)] = 0.50*(GlobalMeshU[GU(i,j,k,0)] + GlobalMeshU[GU(i+1,j,k,0)]);
			}
		}
	}

	// Staggered W mesh
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){
			GlobalMeshW[GW(i,j,NZ,0)] = 0.50*(GlobalMeshU[GU(i,j,NZ-1,0)] + GlobalMeshU[GU(i+1,j,NZ-1,0)]);
			for (k = 0; k < NZ; k++){
				GlobalMeshW[GW(i,j,k,0)] = 0.50*(GlobalMeshU[GU(i,j,k,0)] + GlobalMeshU[GU(i+1,j,k,0)]);
			}
		}
	}


	// Coordinates Y

	if(OptionY == 1){ // Regular

		// Staggered V mesh
		for (i = 0; i < NX; i++){
			for (j = 0; j < NY + 1; j++){
				for (k = 0; k < NZ; k++){
					GlobalMeshV[GV(i,j,k,1)] = j*(Ydominio/ny);
				}
			}
		}

	}
	else if(OptionY == 2){ // Hyperbolic Tangent

		// Staggered V mesh
		for (i = 0; i < NX; i++){
			for (j = 0; j < NY + 1; j++){
				for (k = 0; k < NZ; k++){
					J = j;
					GlobalMeshV[GV(i,j,k,1)] = (Ydominio/2.0)*(1.0 + (tanh(SFY*((2.0*J - ny)/ny)) + tanh(SFY))/tanh(SFY) - 1.0);
				}
			}
		}

	}

	// Collocated Pressure mesh
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){
			for (k = 0; k < NZ; k++){
				GlobalMeshP[GP(i,j,k,1)] = 0.50*(GlobalMeshV[GV(i,j,k,1)] + GlobalMeshV[GV(i,j+1,k,1)]);
			}
		}
	}

	// Staggered U mesh
	for (k = 0; k < NZ; k++){
		for (j = 0; j < NY; j++){
			GlobalMeshU[GU(NX,j,k,1)] = 0.50*(GlobalMeshV[GV(NX-1,j,k,1)] + GlobalMeshV[GV(NX-1,j+1,k,1)]);
			for (i = 0; i < NX; i++){
				GlobalMeshU[GU(i,j,k,1)] = 0.50*(GlobalMeshV[GV(i,j,k,1)] + GlobalMeshV[GV(i,j+1,k,1)]);
			}
		}
	}

	// Staggered W mesh
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){
			GlobalMeshW[GW(i,j,NZ,1)] = 0.50*(GlobalMeshV[GV(i,j,NZ-1,1)] + GlobalMeshV[GV(i,j+1,NZ-1,1)]);
			for (k = 0; k < NZ; k++){
				GlobalMeshW[GW(i,j,k,1)] = 0.50*(GlobalMeshV[GV(i,j,k,1)] + GlobalMeshV[GV(i,j+1,k,1)]);
			}
		}
	}  
 

	// Coordinates Z

	if(OptionZ == 1){ // Regular

		// Staggered W mesh
		for (i = 0; i < NX; i++){
			for (j = 0; j < NY; j++){
				for (k = 0; k < NZ + 1; k++){
					GlobalMeshW[GW(i,j,k,2)] = k*(Zdominio/nz);
				}
			}
		}

	}
	else if(OptionZ == 2){ // Hyperbolic Tangent

		// Staggered W mesh
		for (i = 0; i < NX; i++){
			for (j = 0; j < NY; j++){
				for (k = 0; k < NZ + 1; k++){
					K = k;
					GlobalMeshW[GW(i,j,k,2)] = (Zdominio/2.0)*(1.0 + (tanh(SFZ*((2.0*K - nz)/nz)) + tanh(SFZ))/tanh(SFZ) - 1.0);
				}
			}
		}

	}

	// Collocated Pressure mesh
	for (i = 0; i < NX; i++){
		for (j = 0; j < NY; j++){
			for (k = 0; k < NZ; k++){
				GlobalMeshP[GP(i,j,k,2)] = 0.50*(GlobalMeshW[GW(i,j,k,2)] + GlobalMeshW[GW(i,j,k+1,2)]);
			}
		}
	}

	// Staggered U mesh
	for (k = 0; k < NZ; k++){
		for (j = 0; j < NY; j++){
			GlobalMeshU[GU(NX,j,k,2)] = 0.50*(GlobalMeshW[GW(NX-1,j,k,2)] + GlobalMeshW[GW(NX-1,j,k+1,2)]);
			for (i = 0; i < NX; i++){
				GlobalMeshU[GU(i,j,k,2)] = 0.50*(GlobalMeshW[GW(i,j,k,2)] + GlobalMeshW[GW(i,j,k+1,2)]);
			}
		}
	}

	// Staggered V mesh
	for (i = 0; i < NX; i++){
		for (k = 0; k < NZ; k++){
			GlobalMeshV[GV(i,NY,k,2)] = 0.50*(GlobalMeshW[GW(i,NY-1,k,2)] + GlobalMeshW[GW(i,NY-1,k+1,2)]);
			for (j = 0; j < NY; j++){
				GlobalMeshV[GV(i,j,k,2)] = 0.50*(GlobalMeshW[GW(i,j,k,2)] + GlobalMeshW[GW(i,j,k+1,2)]);
			}
		}
	}

}