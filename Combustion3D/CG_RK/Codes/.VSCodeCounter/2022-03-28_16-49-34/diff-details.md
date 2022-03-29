# Diff Details

Date : 2022-03-28 16:49:34

Directory /home/sergiogus/Documents/Github/TFM/Combustion3D/CG_RK/Codes/CppCodes

Total : 38 files,  5603 codes, 1021 comments, 1603 blanks, all 8227 lines

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
| [CppCodes/Solver.cpp](/CppCodes/Solver.cpp) | C++ | 78 | 16 | 27 | 121 |
| [CppCodes/Solver_CFD_BoundaryConditions_Channel.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Channel.cpp) | C++ | 191 | 26 | 36 | 253 |
| [CppCodes/Solver_CFD_BoundaryConditions_Channel_T.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Channel_T.cpp) | C++ | 105 | 21 | 23 | 149 |
| [CppCodes/Solver_CFD_BoundaryConditions_Differentially.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Differentially.cpp) | C++ | 197 | 26 | 39 | 262 |
| [CppCodes/Solver_CFD_BoundaryConditions_Differentially_T.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Differentially_T.cpp) | C++ | 92 | 19 | 20 | 131 |
| [CppCodes/Solver_CFD_BoundaryConditions_Driven.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Driven.cpp) | C++ | 197 | 26 | 39 | 262 |
| [CppCodes/Solver_CFD_Energy.cpp](/CppCodes/Solver_CFD_Energy.cpp) | C++ | 759 | 64 | 219 | 1,042 |
| [CppCodes/Solver_CFD_Memory.cpp](/CppCodes/Solver_CFD_Memory.cpp) | C++ | 136 | 20 | 72 | 228 |
| [CppCodes/Solver_CFD_MomentumConvection.cpp](/CppCodes/Solver_CFD_MomentumConvection.cpp) | C++ | 614 | 33 | 263 | 910 |
| [CppCodes/Solver_CFD_MomentumDiffusion.cpp](/CppCodes/Solver_CFD_MomentumDiffusion.cpp) | C++ | 365 | 29 | 45 | 439 |
| [CppCodes/Solver_CFD_PoissonCoeffs.cpp](/CppCodes/Solver_CFD_PoissonCoeffs.cpp) | C++ | 167 | 31 | 42 | 240 |
| [CppCodes/Solver_CFD_RK3_Integration.cpp](/CppCodes/Solver_CFD_RK3_Integration.cpp) | C++ | 150 | 44 | 54 | 248 |
| [CppCodes/Solver_CFD_RK4_Integration.cpp](/CppCodes/Solver_CFD_RK4_Integration.cpp) | C++ | 195 | 56 | 69 | 320 |
| [CppCodes/Solver_CFD_Utilities.cpp](/CppCodes/Solver_CFD_Utilities.cpp) | C++ | 332 | 66 | 100 | 498 |
| [CppCodes/Solver_RunSolver.cpp](/CppCodes/Solver_RunSolver.cpp) | C++ | 30 | 293 | 13 | 336 |
| [CppCodes/Solver_Species_BoundaryConditions.cpp](/CppCodes/Solver_Species_BoundaryConditions.cpp) | C++ | 101 | 23 | 24 | 148 |
| [CppCodes/Solver_Species_Chemistry.cpp](/CppCodes/Solver_Species_Chemistry.cpp) | C++ | 0 | 3 | 2 | 5 |
| [CppCodes/Solver_Species_Convection.cpp](/CppCodes/Solver_Species_Convection.cpp) | C++ | 399 | 31 | 153 | 583 |
| [CppCodes/Solver_Species_Diffusion.cpp](/CppCodes/Solver_Species_Diffusion.cpp) | C++ | 371 | 36 | 64 | 471 |
| [CppCodes/Solver_Species_JANAF.cpp](/CppCodes/Solver_Species_JANAF.cpp) | C++ | 92 | 9 | 24 | 125 |
| [CppCodes/Solver_Species_Memory.cpp](/CppCodes/Solver_Species_Memory.cpp) | C++ | 66 | 18 | 33 | 117 |
| [CppCodes/Solver_Species_Read.cpp](/CppCodes/Solver_Species_Read.cpp) | C++ | 80 | 6 | 19 | 105 |
| [CppCodes/Solver_Species_Utilities.cpp](/CppCodes/Solver_Species_Utilities.cpp) | C++ | 79 | 12 | 21 | 112 |
| [HeaderCodes/Memory.h](/HeaderCodes/Memory.h) | C++ | -12 | -6 | -8 | -26 |
| [HeaderCodes/Mesher.h](/HeaderCodes/Mesher.h) | C++ | -68 | -13 | -28 | -109 |
| [HeaderCodes/Parallel.h](/HeaderCodes/Parallel.h) | C++ | -37 | -11 | -17 | -65 |
| [HeaderCodes/PostProcessing.h](/HeaderCodes/PostProcessing.h) | C++ | -47 | -4 | -19 | -70 |
| [HeaderCodes/ReadData.h](/HeaderCodes/ReadData.h) | C++ | -22 | -4 | -12 | -38 |
| [HeaderCodes/Solver.h](/HeaderCodes/Solver.h) | C++ | -250 | -26 | -94 | -370 |

[summary](results.md) / [details](details.md) / [diff summary](diff.md) / diff details