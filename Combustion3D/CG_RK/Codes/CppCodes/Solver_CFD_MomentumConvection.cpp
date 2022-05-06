//------------------------------------------------------------------------------------------------//
//                CPP FILE FOR N-S MOMENTUM EQUATION CONVECTION CALCULATIONS                      //
//------------------------------------------------------------------------------------------------//

// Function to calculate the convective term of Velocity U
void Solver::Get_ConvectionU(Mesher MESH, double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;
double uW_pred, uE_pred;
double uW, uE, vS, vN, uS, uN, wH, wT, uH, uT;

    // Center
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1] - 1; j++){
            for (k = 1; k < NZ - 1; k++){
                
                // Core
                uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)]); 
                uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)]);

                uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U_Matrix[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], EsquemaLargo);
				uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U_Matrix[LU(i+2,j,k,0)], EsquemaLargo);

                vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)]);
                vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V_Matrix[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V_Matrix[LV(i,j+1,k,0)]);

                uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U_Matrix[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], EsquemaLargo);
				uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U_Matrix[LU(i,j+2,k,0)], EsquemaLargo);

                wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)]);
				wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W_Matrix[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W_Matrix[LW(i,j,k+1,0)]);

                uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U_Matrix[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], EsquemaLargo);
				uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U_Matrix[LU(i,j,k+2,0)], EsquemaLargo);

				U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

            }
        }
    }

    // Top
    j = NY - 1;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for (k = 1; k < NZ - 1; k++){
            
            uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)]); 
            uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)]);

            uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U_Matrix[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U_Matrix[LU(i+2,j,k,0)], EsquemaLargo);

            vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)]);
            vN = V.Top[TOP(i,j,k)];

            uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U_Matrix[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], EsquemaLargo);
			uN = U.Top[TOP(i,j,k)];

            wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)]);
			wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W_Matrix[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W_Matrix[LW(i,j,k+1,0)]);

            uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U_Matrix[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], EsquemaLargo);
			uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U_Matrix[LU(i,j,k+2,0)], EsquemaLargo);

			U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

            
        }
    }

    // Bottom 
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0];
        for (k = 1; k < NZ - 1; k++){
                
            uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)]); 
            uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)]);

            uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U_Matrix[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U_Matrix[LU(i+2,j,k,0)], EsquemaLargo);

            vS = V.Bottom[BOTTOM(i,j,k)];
            vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V_Matrix[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V_Matrix[LV(i,j+1,k,0)]);

            uS = U.Bottom[BOTTOM(i,j,k)];
			uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U_Matrix[LU(i,j+2,k,0)], EsquemaLargo);

            wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)]);
			wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W_Matrix[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W_Matrix[LW(i,j,k+1,0)]);

            uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U_Matrix[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], EsquemaLargo);
			uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U_Matrix[LU(i,j,k+2,0)], EsquemaLargo);

			U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

        }
    }

    // Here
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1] - 1; j++){
                
            uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)]); 
            uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)]);

            uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U_Matrix[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U_Matrix[LU(i+2,j,k,0)], EsquemaLargo);

            vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)]);
            vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V_Matrix[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V_Matrix[LV(i,j+1,k,0)]);

            uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U_Matrix[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], EsquemaLargo);
			uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U_Matrix[LU(i,j+2,k,0)], EsquemaLargo);

            wH = W.Here[HERE(i,j,k)];
			wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W_Matrix[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W_Matrix[LW(i,j,k+1,0)]);

            uH = U.Here[HERE(i,j,k)];
			uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U_Matrix[LU(i,j,k+2,0)], EsquemaLargo);

			U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

        }
    }

    // There
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        for(j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1] - 1; j++){
                
            uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)]); 
            uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)]);

            uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U_Matrix[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U_Matrix[LU(i+2,j,k,0)], EsquemaLargo);

            vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)]);
            vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V_Matrix[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V_Matrix[LV(i,j+1,k,0)]);

            uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U_Matrix[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], EsquemaLargo);
			uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U_Matrix[LU(i,j+2,k,0)], EsquemaLargo);

            wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)]);
			wT = W.There[THERE(i,j,k)];

            uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U_Matrix[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], EsquemaLargo);
			uT = U.There[THERE(i,j,k)];

			U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

        }
    }

    // Bottom Here Corner
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0];

        uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)]); 
        uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)]);

        uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U_Matrix[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], EsquemaLargo);
		uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U_Matrix[LU(i+2,j,k,0)], EsquemaLargo);

        vS = V.Bottom[BOTTOM(i,j,k)];
        vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V_Matrix[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V_Matrix[LV(i,j+1,k,0)]);

        uS = U.Bottom[BOTTOM(i,j,k)];
		uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U_Matrix[LU(i,j+2,k,0)], EsquemaLargo);

        wH = W.Here[HERE(i,j,k)];
		wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W_Matrix[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W_Matrix[LW(i,j,k+1,0)]);

        uH = U.Here[HERE(i,j,k)];
		uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U_Matrix[LU(i,j,k+2,0)], EsquemaLargo);

		U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

    }

    // Bottom There Corner
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
        j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0];

        uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)]); 
        uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)]);

        uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U_Matrix[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], EsquemaLargo);
		uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U_Matrix[LU(i+2,j,k,0)], EsquemaLargo);

        vS = V.Bottom[BOTTOM(i,j,k)];
        vN = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j+1,k,0)], V_Matrix[LV(i-1,j+1,k,0)], MESH.MV[LV(i,j+1,k,0)], V_Matrix[LV(i,j+1,k,0)]);

        uS = U.Bottom[BOTTOM(i,j,k)];
		uN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], MESH.MU[LU(i,j+2,k,1)], U_Matrix[LU(i,j+2,k,0)], EsquemaLargo);

        wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)]);
		wT = W.There[THERE(i,j,k)];

        uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U_Matrix[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], EsquemaLargo);
		uT = U.There[THERE(i,j,k)];

		U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

    }

    // Top Here Corner
    j = NY - 1;
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
                
        uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)]); 
        uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)]);

        uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U_Matrix[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], EsquemaLargo);
		uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U_Matrix[LU(i+2,j,k,0)], EsquemaLargo);

        vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)]);
        vN = V.Top[TOP(i,j,k)];

        uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U_Matrix[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], EsquemaLargo);
		uN = U.Top[TOP(i,j,k)];

        wH = W.Here[HERE(i,j,k)];
		wT = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k+1,0)], W_Matrix[LW(i-1,j,k+1,0)], MESH.MW[LW(i,j,k+1,0)], W_Matrix[LW(i,j,k+1,0)]);

        uH = U.Here[HERE(i,j,k)];
		uT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], MESH.MU[LU(i,j,k+2,2)], U_Matrix[LU(i,j,k+2,0)], EsquemaLargo);

		U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

    }

    // Top There Corner
    j = NY - 1;
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango] + 1; i++){
                
            uW_pred = Interpolacion(MESH.MP[LP(i-1,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)]); 
            uE_pred = Interpolacion(MESH.MP[LP(i,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)]);

            uW = ConvectiveScheme(MESH.MP[LP(i-1,j,k,0)], uW_pred, MESH.MU[LU(i-2,j,k,0)], U_Matrix[LU(i-2,j,k,0)], MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], EsquemaLargo);
			uE = ConvectiveScheme(MESH.MP[LP(i,j,k,0)], uE_pred, MESH.MU[LU(i-1,j,k,0)], U_Matrix[LU(i-1,j,k,0)], MESH.MU[LU(i,j,k,0)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i+1,j,k,0)], U_Matrix[LU(i+1,j,k,0)], MESH.MU[LU(i+2,j,k,0)], U_Matrix[LU(i+2,j,k,0)], EsquemaLargo);

            vS = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)]);
            vN = V.Top[TOP(i,j,k)];

            uS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MU[LU(i,j-2,k,1)], U_Matrix[LU(i,j-2,k,0)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j+1,k,1)], U_Matrix[LU(i,j+1,k,0)], EsquemaLargo);
			uN = U.Top[TOP(i,j,k)];

            wH = Interpolacion(MESH.MU[LU(i,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)]);
			wT = W.There[THERE(i,j,k)];

            uH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MU[LU(i,j,k-2,2)], U_Matrix[LU(i,j,k-2,0)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)], MESH.MU[LU(i,j,k+1,2)], U_Matrix[LU(i,j,k+1,0)], EsquemaLargo);
			uT = U.There[THERE(i,j,k)];

			U.Convective[LU(i,j,k,0)] = (1.0 / MESH.VolMU[LU(i,j,k,0)])*(
										  + MESH.SupMU[LU(i,j,k,0)] * (uE * uE - uW * uW)
                                          + MESH.SupMU[LU(i,j,k,1)] * (uN * vN - uS * vS)
                                          + MESH.SupMU[LU(i,j,k,2)] * (uT * wT - uH * wH)
										  );

    }


    if (Rango == 0){
        i = 0;
        for(j = 0; j < NY; j++){
		    for(k = 0; k < NZ; k++){	
				U.Convective[LU(0,j,k,0)] = 0.0;
			}
		}
    }
    else if (Rango == Procesos - 1){
        i = NX;
        for(j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1]; j++){
		    for(k = 0; k < NZ; k++){	
				U.Convective[LU(NX,j,k,0)] = 0.0;
			}
		}
    }

    if (Int_Left){
		i = NX_1;
		for(j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			for(k = 0; k < NZ; k++){	
				U.Convective[LU(NX_1,j,k,0)] = 0.0;
			}
		}
	}

	if (Int_Right){
		i = NX_1 + NX_2;
		for(j = MESH.NY_ColumnMU[i + Halo - Ix[Rango]][0]; j < MESH.NY_ColumnMU[i + Halo - Ix[Rango]][1] - (NY_3 + NY_4); j++){
			for(k = 0; k < NZ; k++){	
				U.Convective[LU(NX_1 + NX_2,j,k,0)] = 0.0;
			}
		}
	}

}

// Function to calculate the convective term of Velocity V
void Solver::Get_ConvectionV(Mesher MESH, double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;
double vS_pred, vN_pred;
double uW, uE, vW, vE, vS, vN, wH, wT, vH, vT;

    // Center
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){
            for (k = 1; k < NZ - 1; k++){

                // Core
                uW = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)]);
				uE = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i+1,j-1,k,1)], U_Matrix[LU(i+1,j-1,k,0)], MESH.MU[LU(i+1,j,k,1)], U_Matrix[LU(i+1,j,k,0)]);

                vW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MV[LV(i-2,j,k,0)], V_Matrix[LV(i-2,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], EsquemaLargo);
				vE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], MESH.MV[LV(i+2,j,k,0)], V_Matrix[LV(i+2,j,k,0)], EsquemaLargo);

				vS_pred = Interpolacion(MESH.MP[LP(i,j-1,k,1)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[LP(i,j,k,1)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)]);

                vS = ConvectiveScheme(MESH.MP[LP(i,j-1,k,1)], vS_pred, MESH.MV[LV(i,j-2,k,1)], V_Matrix[LV(i,j-2,k,0)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[LP(i,j,k,1)], vN_pred, MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], MESH.MV[LV(i,j+2,k,1)], V_Matrix[LV(i,j+2,k,0)], EsquemaLargo);

                wH = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)]);
				wT = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k+1,1)], W_Matrix[LW(i,j-1,k+1,0)], MESH.MW[LW(i,j,k+1,1)], W_Matrix[LW(i,j,k+1,0)]);

                vH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MV[LV(i,j,k-2,2)], V_Matrix[LV(i,j,k-2,0)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], EsquemaLargo);
				vT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], MESH.MV[LV(i,j,k+2,2)], V_Matrix[LV(i,j,k+2,0)], EsquemaLargo);

                // Internal Left
				if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                    uE = 0.0;
                    vE = 0.0;
                }
                // Internal Right
				else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                    uW = 0.0;
                    vW = 0.0;
                }

                V.Convective[LV(i,j,k,0)] = (1.0 / MESH.VolMV[LV(i,j,k,0)])*(
											 + MESH.SupMV[LV(i,j,k,0)] * (vE * uE - vW * uW) 
											 + MESH.SupMV[LV(i,j,k,1)] * (vN * vN - vS * vS)
											 + MESH.SupMV[LV(i,j,k,2)] * (vT * wT - vH * wH)
											 );

            }
        }
    }

    // Here
    k = 0;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){

            // Core
            uW = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)]);
			uE = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i+1,j-1,k,1)], U_Matrix[LU(i+1,j-1,k,0)], MESH.MU[LU(i+1,j,k,1)], U_Matrix[LU(i+1,j,k,0)]);

            vW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MV[LV(i-2,j,k,0)], V_Matrix[LV(i-2,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], EsquemaLargo);
			vE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], MESH.MV[LV(i+2,j,k,0)], V_Matrix[LV(i+2,j,k,0)], EsquemaLargo);

			vS_pred = Interpolacion(MESH.MP[LP(i,j-1,k,1)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)]);
			vN_pred = Interpolacion(MESH.MP[LP(i,j,k,1)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)]);

            vS = ConvectiveScheme(MESH.MP[LP(i,j-1,k,1)], vS_pred, MESH.MV[LV(i,j-2,k,1)], V_Matrix[LV(i,j-2,k,0)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], EsquemaLargo);
			vN = ConvectiveScheme(MESH.MP[LP(i,j,k,1)], vN_pred, MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], MESH.MV[LV(i,j+2,k,1)], V_Matrix[LV(i,j+2,k,0)], EsquemaLargo);

            wH = W.Here[HERE(i,j,k)];
			wT = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k+1,1)], W_Matrix[LW(i,j-1,k+1,0)], MESH.MW[LW(i,j,k+1,1)], W_Matrix[LW(i,j,k+1,0)]);

            vH = V.Here[HERE(i,j,k)];
			vT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], MESH.MV[LV(i,j,k+2,2)], V_Matrix[LV(i,j,k+2,0)], EsquemaLargo);

            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                uE = 0.0;
                vE = 0.0;
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                uW = 0.0;
                vW = 0.0;
            }

            V.Convective[LV(i,j,k,0)] = (1.0 / MESH.VolMV[LV(i,j,k,0)])*(
											 + MESH.SupMV[LV(i,j,k,0)] * (vE * uE - vW * uW) 
											 + MESH.SupMV[LV(i,j,k,1)] * (vN * vN - vS * vS)
											 + MESH.SupMV[LV(i,j,k,2)] * (vT * wT - vH * wH)
											 );

        }
    }

    // There
    k = NZ - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){

            // Core
            uW = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)]);
			uE = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i+1,j-1,k,1)], U_Matrix[LU(i+1,j-1,k,0)], MESH.MU[LU(i+1,j,k,1)], U_Matrix[LU(i+1,j,k,0)]);

            vW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MV[LV(i-2,j,k,0)], V_Matrix[LV(i-2,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], EsquemaLargo);
			vE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], MESH.MV[LV(i+2,j,k,0)], V_Matrix[LV(i+2,j,k,0)], EsquemaLargo);

			vS_pred = Interpolacion(MESH.MP[LP(i,j-1,k,1)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)]);
			vN_pred = Interpolacion(MESH.MP[LP(i,j,k,1)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)]);

            vS = ConvectiveScheme(MESH.MP[LP(i,j-1,k,1)], vS_pred, MESH.MV[LV(i,j-2,k,1)], V_Matrix[LV(i,j-2,k,0)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], EsquemaLargo);
			vN = ConvectiveScheme(MESH.MP[LP(i,j,k,1)], vN_pred, MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], MESH.MV[LV(i,j+2,k,1)], V_Matrix[LV(i,j+2,k,0)], EsquemaLargo);

            wH = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)]);
			wT = W.There[THERE(i,j,k)];

            vH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MV[LV(i,j,k-2,2)], V_Matrix[LV(i,j,k-2,0)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], EsquemaLargo);
			vT = V.There[THERE(i,j,k)];

            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                uE = 0.0;
                vE = 0.0;
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                uW = 0.0;
                vW = 0.0;
            }

            V.Convective[LV(i,j,k,0)] = (1.0 / MESH.VolMV[LV(i,j,k,0)])*(
											 + MESH.SupMV[LV(i,j,k,0)] * (vE * uE - vW * uW) 
											 + MESH.SupMV[LV(i,j,k,1)] * (vN * vN - vS * vS)
											 + MESH.SupMV[LV(i,j,k,2)] * (vT * wT - vH * wH)
											 );

        }
    }

    if (Rango == 0){
        i = 0;

        // Center
        for (j = 1; j < NY; j++){
            for (k = 1; k < NZ - 1; k++){

                uW = U.Left[LEFT(i,j,k)];
				uE = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i+1,j-1,k,1)], U_Matrix[LU(i+1,j-1,k,0)], MESH.MU[LU(i+1,j,k,1)], U_Matrix[LU(i+1,j,k,0)]);

                vW = V.Left[LEFT(i,j,k)];
				vE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], MESH.MV[LV(i+2,j,k,0)], V_Matrix[LV(i+2,j,k,0)], EsquemaLargo);

				vS_pred = Interpolacion(MESH.MP[LP(i,j-1,k,1)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[LP(i,j,k,1)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)]);

                vS = ConvectiveScheme(MESH.MP[LP(i,j-1,k,1)], vS_pred, MESH.MV[LV(i,j-2,k,1)], V_Matrix[LV(i,j-2,k,0)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[LP(i,j,k,1)], vN_pred, MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], MESH.MV[LV(i,j+2,k,1)], V_Matrix[LV(i,j+2,k,0)], EsquemaLargo);

                wH = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)]);
				wT = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k+1,1)], W_Matrix[LW(i,j-1,k+1,0)], MESH.MW[LW(i,j,k+1,1)], W_Matrix[LW(i,j,k+1,0)]);

                vH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MV[LV(i,j,k-2,2)], V_Matrix[LV(i,j,k-2,0)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], EsquemaLargo);
				vT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], MESH.MV[LV(i,j,k+2,2)], V_Matrix[LV(i,j,k+2,0)], EsquemaLargo);

                V.Convective[LV(i,j,k,0)] = (1.0 / MESH.VolMV[LV(i,j,k,0)])*(
											 + MESH.SupMV[LV(i,j,k,0)] * (vE * uE - vW * uW) 
											 + MESH.SupMV[LV(i,j,k,1)] * (vN * vN - vS * vS)
											 + MESH.SupMV[LV(i,j,k,2)] * (vT * wT - vH * wH)
											 );

            }
        }

        // Here
        k = 0;
        for (j = 1; j < NY; j++){

            uW = U.Left[LEFT(i,j,k)];
			uE = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i+1,j-1,k,1)], U_Matrix[LU(i+1,j-1,k,0)], MESH.MU[LU(i+1,j,k,1)], U_Matrix[LU(i+1,j,k,0)]);

            vW = V.Left[LEFT(i,j,k)];
			vE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], MESH.MV[LV(i+2,j,k,0)], V_Matrix[LV(i+2,j,k,0)], EsquemaLargo);

			vS_pred = Interpolacion(MESH.MP[LP(i,j-1,k,1)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)]);
			vN_pred = Interpolacion(MESH.MP[LP(i,j,k,1)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)]);

            vS = ConvectiveScheme(MESH.MP[LP(i,j-1,k,1)], vS_pred, MESH.MV[LV(i,j-2,k,1)], V_Matrix[LV(i,j-2,k,0)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], EsquemaLargo);
			vN = ConvectiveScheme(MESH.MP[LP(i,j,k,1)], vN_pred, MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], MESH.MV[LV(i,j+2,k,1)], V_Matrix[LV(i,j+2,k,0)], EsquemaLargo);

            wH = W.Here[HERE(i,j,k)];
			wT = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k+1,1)], W_Matrix[LW(i,j-1,k+1,0)], MESH.MW[LW(i,j,k+1,1)], W_Matrix[LW(i,j,k+1,0)]);

            vH = V.Here[HERE(i,j,k)];
			vT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], MESH.MV[LV(i,j,k+2,2)], V_Matrix[LV(i,j,k+2,0)], EsquemaLargo);

            V.Convective[LV(i,j,k,0)] = (1.0 / MESH.VolMV[LV(i,j,k,0)])*(
											 + MESH.SupMV[LV(i,j,k,0)] * (vE * uE - vW * uW) 
											 + MESH.SupMV[LV(i,j,k,1)] * (vN * vN - vS * vS)
											 + MESH.SupMV[LV(i,j,k,2)] * (vT * wT - vH * wH)
											 );

        }

        // There
        k = NZ - 1;
        for (j = 1; j < NY; j++){

            uW = U.Left[LEFT(i,j,k)];
			uE = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i+1,j-1,k,1)], U_Matrix[LU(i+1,j-1,k,0)], MESH.MU[LU(i+1,j,k,1)], U_Matrix[LU(i+1,j,k,0)]);

            vW = V.Left[LEFT(i,j,k)];
			vE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], MESH.MV[LV(i+2,j,k,0)], V_Matrix[LV(i+2,j,k,0)], EsquemaLargo);

			vS_pred = Interpolacion(MESH.MP[LP(i,j-1,k,1)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)]);
			vN_pred = Interpolacion(MESH.MP[LP(i,j,k,1)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)]);

            vS = ConvectiveScheme(MESH.MP[LP(i,j-1,k,1)], vS_pred, MESH.MV[LV(i,j-2,k,1)], V_Matrix[LV(i,j-2,k,0)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], EsquemaLargo);
			vN = ConvectiveScheme(MESH.MP[LP(i,j,k,1)], vN_pred, MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], MESH.MV[LV(i,j+2,k,1)], V_Matrix[LV(i,j+2,k,0)], EsquemaLargo);

            wH = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)]);
			wT = W.There[THERE(i,j,k)];

            vH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MV[LV(i,j,k-2,2)], V_Matrix[LV(i,j,k-2,0)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], EsquemaLargo);
			vT = V.There[THERE(i,j,k)];

            V.Convective[LV(i,j,k,0)] = (1.0 / MESH.VolMV[LV(i,j,k,0)])*(
											 + MESH.SupMV[LV(i,j,k,0)] * (vE * uE - vW * uW) 
											 + MESH.SupMV[LV(i,j,k,1)] * (vN * vN - vS * vS)
											 + MESH.SupMV[LV(i,j,k,2)] * (vT * wT - vH * wH)
											 );

        }

    }
    else if (Rango == Procesos - 1){
        i = NX - 1;

        // Center
        for(j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){
            for (k = 1; k < NZ - 1; k++){

                uW = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)]);
				uE = U.Right[RIGHT(i,j,k)];

                vW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MV[LV(i-2,j,k,0)], V_Matrix[LV(i-2,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], EsquemaLargo);
				vE = V.Right[RIGHT(i,j,k)];

				vS_pred = Interpolacion(MESH.MP[LP(i,j-1,k,1)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)]);
				vN_pred = Interpolacion(MESH.MP[LP(i,j,k,1)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)]);

                vS = ConvectiveScheme(MESH.MP[LP(i,j-1,k,1)], vS_pred, MESH.MV[LV(i,j-2,k,1)], V_Matrix[LV(i,j-2,k,0)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], EsquemaLargo);
				vN = ConvectiveScheme(MESH.MP[LP(i,j,k,1)], vN_pred, MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], MESH.MV[LV(i,j+2,k,1)], V_Matrix[LV(i,j+2,k,0)], EsquemaLargo);

                wH = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)]);
				wT = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k+1,1)], W_Matrix[LW(i,j-1,k+1,0)], MESH.MW[LW(i,j,k+1,1)], W_Matrix[LW(i,j,k+1,0)]);

                vH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MV[LV(i,j,k-2,2)], V_Matrix[LV(i,j,k-2,0)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], EsquemaLargo);
				vT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], MESH.MV[LV(i,j,k+2,2)], V_Matrix[LV(i,j,k+2,0)], EsquemaLargo);

                V.Convective[LV(i,j,k,0)] = (1.0 / MESH.VolMV[LV(i,j,k,0)])*(
											 + MESH.SupMV[LV(i,j,k,0)] * (vE * uE - vW * uW) 
											 + MESH.SupMV[LV(i,j,k,1)] * (vN * vN - vS * vS)
											 + MESH.SupMV[LV(i,j,k,2)] * (vT * wT - vH * wH)
											 );

            }
        }

        // Here
        k = 0;
        for(j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){

            uW = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)]);
			uE = U.Right[RIGHT(i,j,k)];

            vW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MV[LV(i-2,j,k,0)], V_Matrix[LV(i-2,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], EsquemaLargo);
			vE = V.Right[RIGHT(i,j,k)];

			vS_pred = Interpolacion(MESH.MP[LP(i,j-1,k,1)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)]);
			vN_pred = Interpolacion(MESH.MP[LP(i,j,k,1)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)]);

            vS = ConvectiveScheme(MESH.MP[LP(i,j-1,k,1)], vS_pred, MESH.MV[LV(i,j-2,k,1)], V_Matrix[LV(i,j-2,k,0)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], EsquemaLargo);
			vN = ConvectiveScheme(MESH.MP[LP(i,j,k,1)], vN_pred, MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], MESH.MV[LV(i,j+2,k,1)], V_Matrix[LV(i,j+2,k,0)], EsquemaLargo);

            wH = W.Here[HERE(i,j,k)];
			wT = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k+1,1)], W_Matrix[LW(i,j-1,k+1,0)], MESH.MW[LW(i,j,k+1,1)], W_Matrix[LW(i,j,k+1,0)]);

            vH = V.Here[HERE(i,j,k)];
			vT = ConvectiveScheme(MESH.MW[LW(i,j,k+1,2)], wT, MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], MESH.MV[LV(i,j,k+2,2)], V_Matrix[LV(i,j,k+2,0)], EsquemaLargo);

            V.Convective[LV(i,j,k,0)] = (1.0 / MESH.VolMV[LV(i,j,k,0)])*(
											 + MESH.SupMV[LV(i,j,k,0)] * (vE * uE - vW * uW) 
											 + MESH.SupMV[LV(i,j,k,1)] * (vN * vN - vS * vS)
											 + MESH.SupMV[LV(i,j,k,2)] * (vT * wT - vH * wH)
											 );

        }

        // There
        k = NZ - 1;
        for(j = MESH.NY_ColumnMV[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMV[i + Halo - Ix[Rango]][1] - 1; j++){

            uW = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MU[LU(i,j-1,k,1)], U_Matrix[LU(i,j-1,k,0)], MESH.MU[LU(i,j,k,1)], U_Matrix[LU(i,j,k,0)]);
			uE = U.Right[RIGHT(i,j,k)];

            vW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MV[LV(i-2,j,k,0)], V_Matrix[LV(i-2,j,k,0)], MESH.MV[LV(i-1,j,k,0)], V_Matrix[LV(i-1,j,k,0)], MESH.MV[LV(i,j,k,0)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i+1,j,k,0)], V_Matrix[LV(i+1,j,k,0)], EsquemaLargo);
			vE = V.Right[RIGHT(i,j,k)];

			vS_pred = Interpolacion(MESH.MP[LP(i,j-1,k,1)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)]);
			vN_pred = Interpolacion(MESH.MP[LP(i,j,k,1)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)]);

            vS = ConvectiveScheme(MESH.MP[LP(i,j-1,k,1)], vS_pred, MESH.MV[LV(i,j-2,k,1)], V_Matrix[LV(i,j-2,k,0)], MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], EsquemaLargo);
			vN = ConvectiveScheme(MESH.MP[LP(i,j,k,1)], vN_pred, MESH.MV[LV(i,j-1,k,1)], V_Matrix[LV(i,j-1,k,0)], MESH.MV[LV(i,j,k,1)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j+1,k,1)], V_Matrix[LV(i,j+1,k,0)], MESH.MV[LV(i,j+2,k,1)], V_Matrix[LV(i,j+2,k,0)], EsquemaLargo);

            wH = Interpolacion(MESH.MV[LV(i,j,k,1)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)]);
			wT = W.There[THERE(i,j,k)];

            vH = ConvectiveScheme(MESH.MW[LW(i,j,k,2)], wH, MESH.MV[LV(i,j,k-2,2)], V_Matrix[LV(i,j,k-2,0)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)], MESH.MV[LV(i,j,k+1,2)], V_Matrix[LV(i,j,k+1,0)], EsquemaLargo);
			vT = V.There[THERE(i,j,k)];

            V.Convective[LV(i,j,k,0)] = (1.0 / MESH.VolMV[LV(i,j,k,0)])*(
											 + MESH.SupMV[LV(i,j,k,0)] * (vE * uE - vW * uW) 
											 + MESH.SupMV[LV(i,j,k,1)] * (vN * vN - vS * vS)
											 + MESH.SupMV[LV(i,j,k,2)] * (vT * wT - vH * wH)
											 );

        }

    }

}
					
// Function to calculate the convective term of Velocity W
void Solver::Get_ConvectionW(Mesher MESH, double *U_Matrix, double *V_Matrix, double *W_Matrix){
int i, j, k;
double wH_pred, wT_pred;
double wH, wT, uW, uE, wW, wE, vS, vN, wS, wN;

    // Center
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for(j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMW[i + Halo - Ix[Rango]][1] - 1; j++){
            for (k = 1; k < NZ; k++){
                
                // Core
				wH_pred = Interpolacion(MESH.MP[LP(i,j,k-1,2)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[LP(i,j,k,2)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)]);

                wH = ConvectiveScheme(MESH.MP[LP(i,j,k-1,2)], wH_pred, MESH.MW[LW(i,j,k-2,2)], W_Matrix[LW(i,j,k-2,0)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[LP(i,j,k,2)], wT_pred, MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], MESH.MW[LW(i,j,k+2,2)], W_Matrix[LW(i,j,k+2,0)], EsquemaLargo);

				uW = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)]);
				uE = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i+1,j,k-1,2)], U_Matrix[LU(i+1,j,k-1,0)], MESH.MU[LU(i+1,j,k,2)], U_Matrix[LU(i+1,j,k,0)]);

                wW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MW[LW(i-2,j,k,0)], W_Matrix[LW(i-2,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], EsquemaLargo);
				wE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], MESH.MW[LW(i+2,j,k,0)], W_Matrix[LW(i+2,j,k,0)], EsquemaLargo);

				vS = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)]);
				vN = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j+1,k-1,2)], V_Matrix[LV(i,j+1,k-1,0)], MESH.MV[LV(i,j+1,k,2)], V_Matrix[LV(i,j+1,k,0)]);
	
				wS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MW[LW(i,j-2,k,1)], W_Matrix[LW(i,j-2,k,0)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], EsquemaLargo);
				wN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], MESH.MW[LW(i,j+2,k,1)], W_Matrix[LW(i,j+2,k,0)], EsquemaLargo);

                // Internal Left
				if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                    uE = 0.0;
                    wE = 0.0;
                }
                // Internal Right
				else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                    uW = 0.0;
                    wW = 0.0;
                }

				W.Convective[LW(i,j,k,0)] = (1 / MESH.VolMW[LW(i,j,k,0)])*(
											 + MESH.SupMW[LW(i,j,k,0)] * (wE * uE - wW * uW)
											 + MESH.SupMW[LW(i,j,k,1)] * (wN * vN - wS * vS)
											 + MESH.SupMW[LW(i,j,k,2)] * (wT * wT - wH * wH)
											 );

            }
        }
    }

    // Bottom
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][0];
        for (k = 1; k < NZ; k++){
                
            // Core
			wH_pred = Interpolacion(MESH.MP[LP(i,j,k-1,2)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)]);
			wT_pred = Interpolacion(MESH.MP[LP(i,j,k,2)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)]);

            wH = ConvectiveScheme(MESH.MP[LP(i,j,k-1,2)], wH_pred, MESH.MW[LW(i,j,k-2,2)], W_Matrix[LW(i,j,k-2,0)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], EsquemaLargo);
			wT = ConvectiveScheme(MESH.MP[LP(i,j,k,2)], wT_pred, MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], MESH.MW[LW(i,j,k+2,2)], W_Matrix[LW(i,j,k+2,0)], EsquemaLargo);

			uW = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)]);
			uE = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i+1,j,k-1,2)], U_Matrix[LU(i+1,j,k-1,0)], MESH.MU[LU(i+1,j,k,2)], U_Matrix[LU(i+1,j,k,0)]);

            wW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MW[LW(i-2,j,k,0)], W_Matrix[LW(i-2,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], EsquemaLargo);
			wE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], MESH.MW[LW(i+2,j,k,0)], W_Matrix[LW(i+2,j,k,0)], EsquemaLargo);

			vS = V.Bottom[BOTTOM(i,j,k)];
			vN = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j+1,k-1,2)], V_Matrix[LV(i,j+1,k-1,0)], MESH.MV[LV(i,j+1,k,2)], V_Matrix[LV(i,j+1,k,0)]);
	
			wS = W.Bottom[BOTTOM(i,j,k)];
			wN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], MESH.MW[LW(i,j+2,k,1)], W_Matrix[LW(i,j+2,k,0)], EsquemaLargo);

            // Internal Left
			if (i == NX_1 - 1 && j < NY - (NY_3 + NY_4)){
                uE = 0.0;
                wE = 0.0;
            }
            // Internal Right
			else if (i == NX_1 + NX_2 && j < NY - (NY_3 + NY_4)){
                uW = 0.0;
                wW = 0.0;
            }

			W.Convective[LW(i,j,k,0)] = (1 / MESH.VolMW[LW(i,j,k,0)])*(
											 + MESH.SupMW[LW(i,j,k,0)] * (wE * uE - wW * uW)
											 + MESH.SupMW[LW(i,j,k,1)] * (wN * vN - wS * vS)
											 + MESH.SupMW[LW(i,j,k,2)] * (wT * wT - wH * wH)
											 );

        }
    }

    // Top
    j = NY - 1;
    for (i = Ix[Rango]; i < Fx[Rango]; i++){
        for (k = 1; k < NZ; k++){
                
			wH_pred = Interpolacion(MESH.MP[LP(i,j,k-1,2)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)]);
			wT_pred = Interpolacion(MESH.MP[LP(i,j,k,2)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)]);

            wH = ConvectiveScheme(MESH.MP[LP(i,j,k-1,2)], wH_pred, MESH.MW[LW(i,j,k-2,2)], W_Matrix[LW(i,j,k-2,0)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], EsquemaLargo);
			wT = ConvectiveScheme(MESH.MP[LP(i,j,k,2)], wT_pred, MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], MESH.MW[LW(i,j,k+2,2)], W_Matrix[LW(i,j,k+2,0)], EsquemaLargo);

			uW = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)]);
			uE = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i+1,j,k-1,2)], U_Matrix[LU(i+1,j,k-1,0)], MESH.MU[LU(i+1,j,k,2)], U_Matrix[LU(i+1,j,k,0)]);

            wW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MW[LW(i-2,j,k,0)], W_Matrix[LW(i-2,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], EsquemaLargo);
			wE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], MESH.MW[LW(i+2,j,k,0)], W_Matrix[LW(i+2,j,k,0)], EsquemaLargo);

			vS = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)]);
			vN = V.Top[TOP(i,j,k)];
	
			wS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MW[LW(i,j-2,k,1)], W_Matrix[LW(i,j-2,k,0)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], EsquemaLargo);
			wN = W.Top[TOP(i,j,k)];

			W.Convective[LW(i,j,k,0)] = (1 / MESH.VolMW[LW(i,j,k,0)])*(
											 + MESH.SupMW[LW(i,j,k,0)] * (wE * uE - wW * uW)
											 + MESH.SupMW[LW(i,j,k,1)] * (wN * vN - wS * vS)
											 + MESH.SupMW[LW(i,j,k,2)] * (wT * wT - wH * wH)
											 );

        }
    }

    if (Rango == 0){

        i = 0;

        // Center
        for (j = 1; j < NY - 1; j++){
            for (k = 1; k < NZ; k++){
                
				wH_pred = Interpolacion(MESH.MP[LP(i,j,k-1,2)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[LP(i,j,k,2)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)]);

                wH = ConvectiveScheme(MESH.MP[LP(i,j,k-1,2)], wH_pred, MESH.MW[LW(i,j,k-2,2)], W_Matrix[LW(i,j,k-2,0)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[LP(i,j,k,2)], wT_pred, MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], MESH.MW[LW(i,j,k+2,2)], W_Matrix[LW(i,j,k+2,0)], EsquemaLargo);

				uW = U.Left[LEFT(i,j,k)];
				uE = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i+1,j,k-1,2)], U_Matrix[LU(i+1,j,k-1,0)], MESH.MU[LU(i+1,j,k,2)], U_Matrix[LU(i+1,j,k,0)]);

                wW = W.Left[LEFT(i,j,k)];
				wE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], MESH.MW[LW(i+2,j,k,0)], W_Matrix[LW(i+2,j,k,0)], EsquemaLargo);

				vS = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)]);
				vN = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j+1,k-1,2)], V_Matrix[LV(i,j+1,k-1,0)], MESH.MV[LV(i,j+1,k,2)], V_Matrix[LV(i,j+1,k,0)]);
	
				wS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MW[LW(i,j-2,k,1)], W_Matrix[LW(i,j-2,k,0)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], EsquemaLargo);
				wN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], MESH.MW[LW(i,j+2,k,1)], W_Matrix[LW(i,j+2,k,0)], EsquemaLargo);

				W.Convective[LW(i,j,k,0)] = (1 / MESH.VolMW[LW(i,j,k,0)])*(
											 + MESH.SupMW[LW(i,j,k,0)] * (wE * uE - wW * uW)
											 + MESH.SupMW[LW(i,j,k,1)] * (wN * vN - wS * vS)
											 + MESH.SupMW[LW(i,j,k,2)] * (wT * wT - wH * wH)
											 );

            }
        }

        // Bottom
        j = 0;
        for (k = 1; k < NZ; k++){
                
			wH_pred = Interpolacion(MESH.MP[LP(i,j,k-1,2)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)]);
			wT_pred = Interpolacion(MESH.MP[LP(i,j,k,2)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)]);

            wH = ConvectiveScheme(MESH.MP[LP(i,j,k-1,2)], wH_pred, MESH.MW[LW(i,j,k-2,2)], W_Matrix[LW(i,j,k-2,0)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], EsquemaLargo);
			wT = ConvectiveScheme(MESH.MP[LP(i,j,k,2)], wT_pred, MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], MESH.MW[LW(i,j,k+2,2)], W_Matrix[LW(i,j,k+2,0)], EsquemaLargo);

			uW = U.Left[LEFT(i,j,k)];
			uE = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i+1,j,k-1,2)], U_Matrix[LU(i+1,j,k-1,0)], MESH.MU[LU(i+1,j,k,2)], U_Matrix[LU(i+1,j,k,0)]);

            wW = W.Left[LEFT(i,j,k)];
			wE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], MESH.MW[LW(i+2,j,k,0)], W_Matrix[LW(i+2,j,k,0)], EsquemaLargo);

			vS = V.Bottom[BOTTOM(i,j,k)];
			vN = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j+1,k-1,2)], V_Matrix[LV(i,j+1,k-1,0)], MESH.MV[LV(i,j+1,k,2)], V_Matrix[LV(i,j+1,k,0)]);
	
			wS = W.Bottom[BOTTOM(i,j,k)];
			wN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], MESH.MW[LW(i,j+2,k,1)], W_Matrix[LW(i,j+2,k,0)], EsquemaLargo);

			W.Convective[LW(i,j,k,0)] = (1 / MESH.VolMW[LW(i,j,k,0)])*(
											 + MESH.SupMW[LW(i,j,k,0)] * (wE * uE - wW * uW)
											 + MESH.SupMW[LW(i,j,k,1)] * (wN * vN - wS * vS)
											 + MESH.SupMW[LW(i,j,k,2)] * (wT * wT - wH * wH)
											 );

        }

        // Top
        j = NY - 1;
        for (k = 1; k < NZ; k++){
                
			wH_pred = Interpolacion(MESH.MP[LP(i,j,k-1,2)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)]);
			wT_pred = Interpolacion(MESH.MP[LP(i,j,k,2)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)]);

            wH = ConvectiveScheme(MESH.MP[LP(i,j,k-1,2)], wH_pred, MESH.MW[LW(i,j,k-2,2)], W_Matrix[LW(i,j,k-2,0)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], EsquemaLargo);
			wT = ConvectiveScheme(MESH.MP[LP(i,j,k,2)], wT_pred, MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], MESH.MW[LW(i,j,k+2,2)], W_Matrix[LW(i,j,k+2,0)], EsquemaLargo);

			uW = U.Left[LEFT(i,j,k)];
			uE = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i+1,j,k-1,2)], U_Matrix[LU(i+1,j,k-1,0)], MESH.MU[LU(i+1,j,k,2)], U_Matrix[LU(i+1,j,k,0)]);

            wW = W.Left[LEFT(i,j,k)];
			wE = ConvectiveScheme(MESH.MU[LU(i+1,j,k,0)], uE, MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], MESH.MW[LW(i+2,j,k,0)], W_Matrix[LW(i+2,j,k,0)], EsquemaLargo);

			vS = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)]);
			vN = V.Top[TOP(i,j,k)];
	
			wS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MW[LW(i,j-2,k,1)], W_Matrix[LW(i,j-2,k,0)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], EsquemaLargo);
			wN = W.Top[TOP(i,j,k)];

			W.Convective[LW(i,j,k,0)] = (1 / MESH.VolMW[LW(i,j,k,0)])*(
											 + MESH.SupMW[LW(i,j,k,0)] * (wE * uE - wW * uW)
											 + MESH.SupMW[LW(i,j,k,1)] * (wN * vN - wS * vS)
											 + MESH.SupMW[LW(i,j,k,2)] * (wT * wT - wH * wH)
											 );

        }

    }
    else if (Rango == Procesos - 1){

        i = NX - 1;

        // Center
        for(j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][0] + 1; j < MESH.NY_ColumnMW[i + Halo - Ix[Rango]][1] - 1; j++){
            for (k = 1; k < NZ; k++){
                
				wH_pred = Interpolacion(MESH.MP[LP(i,j,k-1,2)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)]);
				wT_pred = Interpolacion(MESH.MP[LP(i,j,k,2)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)]);

                wH = ConvectiveScheme(MESH.MP[LP(i,j,k-1,2)], wH_pred, MESH.MW[LW(i,j,k-2,2)], W_Matrix[LW(i,j,k-2,0)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], EsquemaLargo);
				wT = ConvectiveScheme(MESH.MP[LP(i,j,k,2)], wT_pred, MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], MESH.MW[LW(i,j,k+2,2)], W_Matrix[LW(i,j,k+2,0)], EsquemaLargo);

				uW = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)]);
				uE = U.Right[RIGHT(i,j,k)];

                wW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MW[LW(i-2,j,k,0)], W_Matrix[LW(i-2,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], EsquemaLargo);
				wE = W.Right[RIGHT(i,j,k)];

				vS = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)]);
				vN = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j+1,k-1,2)], V_Matrix[LV(i,j+1,k-1,0)], MESH.MV[LV(i,j+1,k,2)], V_Matrix[LV(i,j+1,k,0)]);
	
				wS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MW[LW(i,j-2,k,1)], W_Matrix[LW(i,j-2,k,0)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], EsquemaLargo);
				wN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], MESH.MW[LW(i,j+2,k,1)], W_Matrix[LW(i,j+2,k,0)], EsquemaLargo);

				W.Convective[LW(i,j,k,0)] = (1 / MESH.VolMW[LW(i,j,k,0)])*(
											 + MESH.SupMW[LW(i,j,k,0)] * (wE * uE - wW * uW)
											 + MESH.SupMW[LW(i,j,k,1)] * (wN * vN - wS * vS)
											 + MESH.SupMW[LW(i,j,k,2)] * (wT * wT - wH * wH)
											 );

            }
        }

        // Bottom
        j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][0];
        for (k = 1; k < NZ; k++){
            
			wH_pred = Interpolacion(MESH.MP[LP(i,j,k-1,2)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)]);
			wT_pred = Interpolacion(MESH.MP[LP(i,j,k,2)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)]);

            wH = ConvectiveScheme(MESH.MP[LP(i,j,k-1,2)], wH_pred, MESH.MW[LW(i,j,k-2,2)], W_Matrix[LW(i,j,k-2,0)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], EsquemaLargo);
			wT = ConvectiveScheme(MESH.MP[LP(i,j,k,2)], wT_pred, MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], MESH.MW[LW(i,j,k+2,2)], W_Matrix[LW(i,j,k+2,0)], EsquemaLargo);

			uW = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)]);
			uE = U.Right[RIGHT(i,j,k)];

            wW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MW[LW(i-2,j,k,0)], W_Matrix[LW(i-2,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], EsquemaLargo);
			wE = W.Right[RIGHT(i,j,k)];

			vS = V.Bottom[BOTTOM(i,j,k)];
			vN = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j+1,k-1,2)], V_Matrix[LV(i,j+1,k-1,0)], MESH.MV[LV(i,j+1,k,2)], V_Matrix[LV(i,j+1,k,0)]);
	
			wS = W.Bottom[BOTTOM(i,j,k)];
			wN = ConvectiveScheme(MESH.MV[LV(i,j+1,k,1)], vN, MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], MESH.MW[LW(i,j+2,k,1)], W_Matrix[LW(i,j+2,k,0)], EsquemaLargo);

			W.Convective[LW(i,j,k,0)] = (1 / MESH.VolMW[LW(i,j,k,0)])*(
											 + MESH.SupMW[LW(i,j,k,0)] * (wE * uE - wW * uW)
											 + MESH.SupMW[LW(i,j,k,1)] * (wN * vN - wS * vS)
											 + MESH.SupMW[LW(i,j,k,2)] * (wT * wT - wH * wH)
											 );

        }

        // Top
        j = MESH.NY_ColumnMW[i + Halo - Ix[Rango]][1] - 1;
        for (k = 1; k < NZ; k++){
                
			wH_pred = Interpolacion(MESH.MP[LP(i,j,k-1,2)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)]);
			wT_pred = Interpolacion(MESH.MP[LP(i,j,k,2)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)]);

            wH = ConvectiveScheme(MESH.MP[LP(i,j,k-1,2)], wH_pred, MESH.MW[LW(i,j,k-2,2)], W_Matrix[LW(i,j,k-2,0)], MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], EsquemaLargo);
			wT = ConvectiveScheme(MESH.MP[LP(i,j,k,2)], wT_pred, MESH.MW[LW(i,j,k-1,2)], W_Matrix[LW(i,j,k-1,0)], MESH.MW[LW(i,j,k,2)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j,k+1,2)], W_Matrix[LW(i,j,k+1,0)], MESH.MW[LW(i,j,k+2,2)], W_Matrix[LW(i,j,k+2,0)], EsquemaLargo);

			uW = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MU[LU(i,j,k-1,2)], U_Matrix[LU(i,j,k-1,0)], MESH.MU[LU(i,j,k,2)], U_Matrix[LU(i,j,k,0)]);
			uE = U.Right[RIGHT(i,j,k)];

            wW = ConvectiveScheme(MESH.MU[LU(i,j,k,0)], uW, MESH.MW[LW(i-2,j,k,0)], W_Matrix[LW(i-2,j,k,0)], MESH.MW[LW(i-1,j,k,0)], W_Matrix[LW(i-1,j,k,0)], MESH.MW[LW(i,j,k,0)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i+1,j,k,0)], W_Matrix[LW(i+1,j,k,0)], EsquemaLargo);
			wE = W.Right[RIGHT(i,j,k)];

			vS = Interpolacion(MESH.MW[LW(i,j,k,2)], MESH.MV[LV(i,j,k-1,2)], V_Matrix[LV(i,j,k-1,0)], MESH.MV[LV(i,j,k,2)], V_Matrix[LV(i,j,k,0)]);
			vN = V.Top[TOP(i,j,k)];
	
			wS = ConvectiveScheme(MESH.MV[LV(i,j,k,1)], vS, MESH.MW[LW(i,j-2,k,1)], W_Matrix[LW(i,j-2,k,0)], MESH.MW[LW(i,j-1,k,1)], W_Matrix[LW(i,j-1,k,0)], MESH.MW[LW(i,j,k,1)], W_Matrix[LW(i,j,k,0)], MESH.MW[LW(i,j+1,k,1)], W_Matrix[LW(i,j+1,k,0)], EsquemaLargo);
			wN = W.Top[TOP(i,j,k)];

			W.Convective[LW(i,j,k,0)] = (1 / MESH.VolMW[LW(i,j,k,0)])*(
											 + MESH.SupMW[LW(i,j,k,0)] * (wE * uE - wW * uW)
											 + MESH.SupMW[LW(i,j,k,1)] * (wN * vN - wS * vS)
											 + MESH.SupMW[LW(i,j,k,2)] * (wT * wT - wH * wH)
											 );

        }

    }

}
				
					