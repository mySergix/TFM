# Diff Details

Date : 2022-03-28 09:20:09

Directory /home/sergiogus/Documents/Github/TFM/Combustion3D/CG_RK/Codes/CppCodes

Total : 53 files,  1687 codes, 205 comments, 446 blanks, all 2338 lines

[summary](results.md) / [details](details.md) / [diff summary](diff.md) / diff details

## Files
| filename | language | code | comment | blank | total |
| :--- | :--- | ---: | ---: | ---: | ---: |
| [CppCodes/Matrix_Index.cpp](/CppCodes/Matrix_Index.cpp) | C++ | 19 | 25 | 24 | 68 |
| [CppCodes/Memory.cpp](/CppCodes/Memory.cpp) | C++ | 14 | 5 | 15 | 34 |
| [CppCodes/Mesher.cpp](/CppCodes/Mesher.cpp) | C++ | 274 | 31 | 68 | 373 |
| [CppCodes/Mesher_Memory.cpp](/CppCodes/Mesher_Memory.cpp) | C++ | 56 | 18 | 22 | 96 |
| [CppCodes/Mesher_NodalCoordinates.cpp](/CppCodes/Mesher_NodalCoordinates.cpp) | C++ | 266 | 41 | 61 | 368 |
| [CppCodes/Parallel.cpp](/CppCodes/Parallel.cpp) | C++ | 74 | 12 | 24 | 110 |
| [CppCodes/Parallel_Communication.cpp](/CppCodes/Parallel_Communication.cpp) | C++ | 162 | 12 | 64 | 238 |
| [CppCodes/PostProcessing.cpp](/CppCodes/PostProcessing.cpp) | C++ | 314 | 25 | 79 | 418 |
| [CppCodes/ReadData.cpp](/CppCodes/ReadData.cpp) | C++ | 64 | 8 | 23 | 95 |
| [CppCodes/Solver.cpp](/CppCodes/Solver.cpp) | C++ | 73 | 15 | 25 | 113 |
| [CppCodes/Solver_CFD_BoundaryConditions_Channel.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Channel.cpp) | C++ | 191 | 26 | 36 | 253 |
| [CppCodes/Solver_CFD_BoundaryConditions_Channel_T.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Channel_T.cpp) | C++ | 105 | 21 | 23 | 149 |
| [CppCodes/Solver_CFD_BoundaryConditions_Differentially.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Differentially.cpp) | C++ | 197 | 26 | 39 | 262 |
| [CppCodes/Solver_CFD_BoundaryConditions_Differentially_T.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Differentially_T.cpp) | C++ | 92 | 19 | 20 | 131 |
| [CppCodes/Solver_CFD_BoundaryConditions_Driven.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Driven.cpp) | C++ | 197 | 26 | 39 | 262 |
| [CppCodes/Solver_CFD_Energy.cpp](/CppCodes/Solver_CFD_Energy.cpp) | C++ | 743 | 62 | 213 | 1,018 |
| [CppCodes/Solver_CFD_Memory.cpp](/CppCodes/Solver_CFD_Memory.cpp) | C++ | 135 | 20 | 72 | 227 |
| [CppCodes/Solver_CFD_MomentumConvection.cpp](/CppCodes/Solver_CFD_MomentumConvection.cpp) | C++ | 614 | 33 | 263 | 910 |
| [CppCodes/Solver_CFD_MomentumDiffusion.cpp](/CppCodes/Solver_CFD_MomentumDiffusion.cpp) | C++ | 365 | 29 | 45 | 439 |
| [CppCodes/Solver_CFD_PoissonCoeffs.cpp](/CppCodes/Solver_CFD_PoissonCoeffs.cpp) | C++ | 167 | 31 | 42 | 240 |
| [CppCodes/Solver_CFD_RK3_Integration.cpp](/CppCodes/Solver_CFD_RK3_Integration.cpp) | C++ | 150 | 44 | 54 | 248 |
| [CppCodes/Solver_CFD_RK4_Integration.cpp](/CppCodes/Solver_CFD_RK4_Integration.cpp) | C++ | 195 | 56 | 69 | 320 |
| [CppCodes/Solver_CFD_Utilities.cpp](/CppCodes/Solver_CFD_Utilities.cpp) | C++ | 332 | 66 | 100 | 498 |
| [CppCodes/Solver_RunSolver.cpp](/CppCodes/Solver_RunSolver.cpp) | C++ | 207 | 43 | 81 | 331 |
| [CppCodes/Solver_Species_BoundaryConditions.cpp](/CppCodes/Solver_Species_BoundaryConditions.cpp) | C++ | 101 | 23 | 24 | 148 |
| [CppCodes/Solver_Species_Convection.cpp](/CppCodes/Solver_Species_Convection.cpp) | C++ | 399 | 31 | 153 | 583 |
| [CppCodes/Solver_Species_Diffusion.cpp](/CppCodes/Solver_Species_Diffusion.cpp) | C++ | 329 | 32 | 47 | 408 |
| [CppCodes/Solver_Species_JANAF.cpp](/CppCodes/Solver_Species_JANAF.cpp) | C++ | 82 | 8 | 20 | 110 |
| [CppCodes/Solver_Species_Memory.cpp](/CppCodes/Solver_Species_Memory.cpp) | C++ | 60 | 16 | 30 | 106 |
| [CppCodes/Solver_Species_Read.cpp](/CppCodes/Solver_Species_Read.cpp) | C++ | 80 | 6 | 19 | 105 |
| [CppCodes/Solver_Species_Utilities.cpp](/CppCodes/Solver_Species_Utilities.cpp) | C++ | 75 | 11 | 21 | 107 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Matrix_Index.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Matrix_Index.cpp) | C++ | -18 | -24 | -23 | -65 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Memory.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Memory.cpp) | C++ | -14 | -5 | -15 | -34 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Mesher.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Mesher.cpp) | C++ | -214 | -26 | -57 | -297 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Mesher_Memory.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Mesher_Memory.cpp) | C++ | -52 | -18 | -21 | -91 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Mesher_NodalCoordinates.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Mesher_NodalCoordinates.cpp) | C++ | -266 | -41 | -61 | -368 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Parallel.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Parallel.cpp) | C++ | -74 | -12 | -24 | -110 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Parallel_Communication.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Parallel_Communication.cpp) | C++ | -144 | -11 | -57 | -212 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/PostProcessing.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/PostProcessing.cpp) | C++ | -210 | -12 | -54 | -276 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/ReadData.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/ReadData.cpp) | C++ | -64 | -8 | -23 | -95 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver.cpp) | C++ | -60 | -10 | -19 | -89 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_BoundaryConditions_Differentially.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_BoundaryConditions_Differentially.cpp) | C++ | -197 | -26 | -39 | -262 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_BoundaryConditions_Differentially_T.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_BoundaryConditions_Differentially_T.cpp) | C++ | -92 | -19 | -20 | -131 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_BoundaryConditions_Driven.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_BoundaryConditions_Driven.cpp) | C++ | -197 | -26 | -39 | -262 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_Energy.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_Energy.cpp) | C++ | -743 | -62 | -213 | -1,018 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_Memory.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_Memory.cpp) | C++ | -132 | -20 | -69 | -221 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_MomentumConvection.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_MomentumConvection.cpp) | C++ | -614 | -33 | -263 | -910 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_MomentumDiffusion.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_MomentumDiffusion.cpp) | C++ | -365 | -29 | -45 | -439 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_PoissonCoeffs.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_PoissonCoeffs.cpp) | C++ | -167 | -31 | -42 | -240 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_RK3_Integration.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_RK3_Integration.cpp) | C++ | -147 | -44 | -54 | -245 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_RK4_Integration.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_RK4_Integration.cpp) | C++ | -191 | -56 | -69 | -316 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_Utilities.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_CFD_Utilities.cpp) | C++ | -294 | -62 | -89 | -445 |
| [/home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_RunSolver.cpp](//home/sergiogus/Documents/Github/TFM/NavierStokes3D/CG_RK3/Codes/CppCodes/Solver_RunSolver.cpp) | C++ | -190 | -41 | -73 | -304 |

[summary](results.md) / [details](details.md) / [diff summary](diff.md) / diff details