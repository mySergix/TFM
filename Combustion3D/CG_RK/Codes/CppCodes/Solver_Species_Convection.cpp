//------------------------------------------------------------------------------------------------//
//                   CPP FILE FOR SPECIES EQUATION CONVECTION CALCULATIONS                        //
//------------------------------------------------------------------------------------------------//

// Function to calculate the convection of each species
void Solver::Get_Species_Convection(Mesher MESH){
int i, j, k, n;
double Yw, Ye, Ys, Yn, Yh, Yt;

    for (n = 0; n < N_Species - 1; n++){

        // Center
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){
                for (k = 1; k < NZ - 1; k++){

                    // Core
                    Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                    Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                    Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                    Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                    Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                    Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                    // Internal Left
				    if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                        Ye = Species[n].Y_Pres[LP(i,j,k,0)];
                    }
                    // Internal Right
				    else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                        Yw = Species[n].Y_Pres[LP(i,j,k,0)];
                    }

                    Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

                }
            }
        }

        // Bottom
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
            for (k = 1; k < NZ - 1; k++){

                // Core
                Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Ys = Species[n].Bottom[BOTTOM(i,j,k)];
                Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                // Internal Left
				if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                    Ye = Species[n].Y_Pres[LP(i,j,k,0)];
                }
                // Internal Right
				else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                    Yw = Species[n].Y_Pres[LP(i,j,k,0)];
                }

                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

            }
        }

        // Top
        j = NY - 1;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for (k = 1; k < NZ - 1; k++){
                Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Yn = Species[n].Top[TOP(i,j,k)];

                Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Top[TOP(i,j,k)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

            }
        }

        // Here
        k = 0;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){

                // Core
                Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Yh = Species[n].Here[HERE(i,j,k)];
                Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                // Internal Left
				if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                    Ye = Species[n].Y_Pres[LP(i,j,k,0)];
                }
                // Internal Right
				else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                    Yw = Species[n].Y_Pres[LP(i,j,k,0)];
                }

                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Here[HERE(i,j,k)])
											 );

            }
        }

        // There
        k = NZ - 1;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){

                // Core
                Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Yt = Species[n].There[THERE(i,j,k)];

                // Internal Left
				if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                    Ye = Species[n].Y_Pres[LP(i,j,k,0)];
                }
                // Internal Right
				else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                    Yw = Species[n].Y_Pres[LP(i,j,k,0)];
                }

                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

            }
        }

        // Bottom Here Corner
        k = 0;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];

            Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ys = Species[n].Bottom[BOTTOM(i,j,k)];
            Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Yh = Species[n].Here[HERE(i,j,k)];
            Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Here[HERE(i,j,k)])
											 );

        }

        // Bottom There Corner
        k = NZ - 1;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];

            Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ys = Species[n].Bottom[BOTTOM(i,j,k)];
            Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Yt = Species[n].There[THERE(i,j,k)];

            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

        }

        // Top Here Corner
        j = NY - 1;
        k = 0;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Yn = Species[n].Top[TOP(i,j,k)];
            
            Yh = Species[n].Here[HERE(i,j,k)];
            Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Here[HERE(i,j,k)])
											 );

        }

        // Top There Corner
        j = NY - 1;
        k = NZ - 1;
        for (i = Ix[Rango]; i < Fx[Rango]; i++){
            Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Yn = Species[n].Top[TOP(i,j,k)];

            Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Yt = Species[n].There[THERE(i,j,k)];

            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

        }
    
        if (Rango == 0){
            i = 0;

            // Center
            for (j = 1; j < NY - 1; j++){
                for (k = 1; k < NZ - 1; k++){
                    Yw = Species[n].Left[LEFT(i,j,k)];
                    Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                    Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                    Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                    Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                    Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                    Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

                }
            }

            // Bottom
            j = 0;
            for (k = 1; k < NZ - 1; k++){
                Yw = Species[n].Left[LEFT(i,j,k)];
                Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Ys = Species[n].Bottom[BOTTOM(i,j,k)];
                Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

            }

            // Top
            j = NY - 1;
            for (k = 1; k < NZ - 1; k++){
                Yw = Species[n].Left[LEFT(i,j,k)];
                Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Yn = Species[n].Top[TOP(i,j,k)];

                Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

            }

            // Here
            k = 0;
            for (j = 1; j < NY - 1; j++){
                Yw = Species[n].Left[LEFT(i,j,k)];
                Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Yh = Species[n].Here[HERE(i,j,k)];
                Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

            }

            // There
            k = NZ - 1;
            for (j = 1; j < NY - 1; j++){
                Yw = Species[n].Left[LEFT(i,j,k)];
                Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
                Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Yt = Species[n].There[THERE(i,j,k)];

                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

            }

            // Bottom Here Corner
            j = 0;
            k = 0;

            Yw = Species[n].Left[LEFT(i,j,k)];
            Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ys = Species[n].Bottom[BOTTOM(i,j,k)];
            Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Yh = Species[n].Here[HERE(i,j,k)];
            Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Here[HERE(i,j,k)])
											 );


            // Bottom There Corner
            j = 0;
            k = NZ - 1;

            Yw = Species[n].Left[LEFT(i,j,k)];
            Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ys = Species[n].Bottom[BOTTOM(i,j,k)];
            Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Yt = Species[n].There[THERE(i,j,k)];

            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

    

            // Top Here Corner
            j = NY - 1;
            k = 0;

            Yw = Species[n].Left[LEFT(i,j,k)];
            Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Yn = Species[n].Top[TOP(i,j,k)];
            
            Yh = Species[n].Here[HERE(i,j,k)];
            Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Here[HERE(i,j,k)])
											 );

    

            // Top There Corner
            j = NY - 1;
            k = NZ - 1;

            Yw = Species[n].Left[LEFT(i,j,k)];
            Ye = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], U.Pres[LU(i+1,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], MESH.MP[LP(i+2,j,k,0)], Species[n].Y_Pres[LP(i+2,j,k,0)], EsquemaLargo);
            
            Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Yn = Species[n].Top[TOP(i,j,k)];
            
            Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Yt = Species[n].There[THERE(i,j,k)];

            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

        }
        else if (Rango == Procesos - 1){
            i = NX - 1;

            // Center
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){
                for (k = 1; k < NZ - 1; k++){
                    Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                    Ye = Species[n].Right[RIGHT(i,j,k)];
                
                    Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                    Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                    Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                    Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                    Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

                }
            }

            // Bottom
            j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
            for (k = 1; k < NZ - 1; k++){
                Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Ye = Species[n].Right[RIGHT(i,j,k)];

                Ys = Species[n].Bottom[BOTTOM(i,j,k)];
                Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);

                Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

            }

            // Top
            j = NY - 1;
            for (k = 1; k < NZ - 1; k++){
                Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Ye = Species[n].Right[RIGHT(i,j,k)];

                Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Yn = Species[n].Top[TOP(i,j,k)];

                Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

            }

            // Here
            k = 0;
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){
                Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Ye = Species[n].Right[RIGHT(i,j,k)];

                Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Yh = Species[n].Here[HERE(i,j,k)];
                Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

            }

            // There
            k = NZ - 1;
            for(j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMP[i + HP - Ix[Rango]][1] - 1; j++){
                Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
                Ye = Species[n].Right[RIGHT(i,j,k)];

                Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
                Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
                Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
                Yt = Species[n].There[THERE(i,j,k)];

                Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

            }

            // Bottom Here Corner
            j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
            k = 0;
    
            Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ye = Species[n].Right[RIGHT(i,j,k)];

            Ys = Species[n].Bottom[BOTTOM(i,j,k)];
            Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Yh = Species[n].Here[HERE(i,j,k)];
            Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Here[HERE(i,j,k)])
											 );

    

            // Bottom There Corner
            j = MESH.NY_ColumnMP[i + HP - Ix[Rango]][0];
            k = NZ - 1;
    
            Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ye = Species[n].Right[RIGHT(i,j,k)];

            Ys = Species[n].Bottom[BOTTOM(i,j,k)];
            Yn = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], V.Pres[LV(i,j+1,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], MESH.MP[LP(i,j+2,k,1)], Species[n].Y_Pres[LP(i,j+2,k,0)], EsquemaLargo);
            
            Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Yt = Species[n].There[THERE(i,j,k)];

            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

    

            // Top Here Corner
            j = NY - 1;
            k = 0;
    
            Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ye = Species[n].Right[RIGHT(i,j,k)];

            Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Yn = Species[n].Top[TOP(i,j,k)];
            
            Yh = Species[n].Here[HERE(i,j,k)];
            Yt = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], W.Pres[LW(i,j,k+1,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], MESH.MP[LP(i,j,k+2,2)], Species[n].Y_Pres[LP(i,j,k+2,0)], EsquemaLargo);
            
            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Here[HERE(i,j,k)])
											 );

    

            // Top There Corner
            j = NY - 1;
            k = NZ - 1;
    
            Yw = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], U.Pres[LU(i,j,k,0)], MESH.MP[LP(i-2,j,k,0)], Species[n].Y_Pres[LP(i-2,j,k,0)], MESH.MP[LP(i-1,j,k,0)], Species[n].Y_Pres[LP(i-1,j,k,0)], MESH.MP[LP(i,j,k,0)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i+1,j,k,0)], Species[n].Y_Pres[LP(i+1,j,k,0)], EsquemaLargo);
            Ye = Species[n].Right[RIGHT(i,j,k)];
        
            Ys = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], V.Pres[LV(i,j,k,0)], MESH.MP[LP(i,j-2,k,1)], Species[n].Y_Pres[LP(i,j-2,k,0)], MESH.MP[LP(i,j-1,k,1)], Species[n].Y_Pres[LP(i,j-1,k,0)], MESH.MP[LP(i,j,k,1)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j+1,k,1)], Species[n].Y_Pres[LP(i,j+1,k,0)], EsquemaLargo);
            Yn = Species[n].Top[TOP(i,j,k)];

            Yh = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], W.Pres[LW(i,j,k,0)], MESH.MP[LP(i,j,k-2,2)], Species[n].Y_Pres[LP(i,j,k-2,0)], MESH.MP[LP(i,j,k-1,2)], Species[n].Y_Pres[LP(i,j,k-1,0)], MESH.MP[LP(i,j,k,2)], Species[n].Y_Pres[LP(i,j,k,0)], MESH.MP[LP(i,j,k+1,2)], Species[n].Y_Pres[LP(i,j,k+1,0)], EsquemaLargo);
            Yt = Species[n].There[THERE(i,j,k)];

            Species[n].Convective[LP(i,j,k,0)] = (1.0 / MESH.VolMP[LP(i,j,k,0)])*(
											 + MESH.SupMP[LP(i,j,k,0)] * (U.Pres[LU(i+1,j,k,0)] * Ye - Yw * U.Pres[LU(i,j,k,0)]) 
											 + MESH.SupMP[LP(i,j,k,1)] * (V.Pres[LV(i,j+1,k,0)] * Yn - Ys * V.Pres[LV(i,j,k,0)])
											 + MESH.SupMP[LP(i,j,k,2)] * (W.Pres[LW(i,j,k+1,0)] * Yt - Yh * W.Pres[LW(i,j,k,0)])
											 );

        }
    }

}

