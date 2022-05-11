# Diff Details

Date : 2022-05-11 10:12:58

Directory /home/sergiogus/Documents/Github/TFM/Combustion3D/CG_RK/Codes/CppCodes

Total : 46 files,  1626 codes, 376 comments, 363 blanks, all 2365 lines

[Summary](results.md) / [Details](details.md) / [Diff Summary](diff.md) / Diff Details

## Files
| filename | language | code | comment | blank | total |
| :--- | :--- | ---: | ---: | ---: | ---: |
| [CppCodes/Matrix_Index.cpp](/CppCodes/Matrix_Index.cpp) | C++ | 20 | 26 | 24 | 70 |
| [CppCodes/Memory.cpp](/CppCodes/Memory.cpp) | C++ | 23 | 6 | 19 | 48 |
| [CppCodes/Mesher.cpp](/CppCodes/Mesher.cpp) | C++ | 479 | 50 | 152 | 681 |
| [CppCodes/Mesher_Memory.cpp](/CppCodes/Mesher_Memory.cpp) | C++ | 66 | 20 | 24 | 110 |
| [CppCodes/Mesher_NodalCoordinates.cpp](/CppCodes/Mesher_NodalCoordinates.cpp) | C++ | 466 | 121 | 86 | 673 |
| [CppCodes/Mesher_Paraview.cpp](/CppCodes/Mesher_Paraview.cpp) | C++ | 37 | 4 | 10 | 51 |
| [CppCodes/Parallel.cpp](/CppCodes/Parallel.cpp) | C++ | 100 | 16 | 32 | 148 |
| [CppCodes/Parallel_Communication.cpp](/CppCodes/Parallel_Communication.cpp) | C++ | 144 | 11 | 57 | 212 |
| [CppCodes/PostProcessing.cpp](/CppCodes/PostProcessing.cpp) | C++ | 163 | 23 | 51 | 237 |
| [CppCodes/ReadData.cpp](/CppCodes/ReadData.cpp) | C++ | 64 | 8 | 23 | 95 |
| [CppCodes/Solver.cpp](/CppCodes/Solver.cpp) | C++ | 77 | 13 | 29 | 119 |
| [CppCodes/Solver_CFD_BoundaryConditions_NonPremixed.cpp](/CppCodes/Solver_CFD_BoundaryConditions_NonPremixed.cpp) | C++ | 268 | 32 | 52 | 352 |
| [CppCodes/Solver_CFD_BoundaryConditions_NonPremixed_T.cpp](/CppCodes/Solver_CFD_BoundaryConditions_NonPremixed_T.cpp) | C++ | 135 | 23 | 25 | 183 |
| [CppCodes/Solver_CFD_BoundaryConditions_Premixed.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Premixed.cpp) | C++ | 277 | 32 | 56 | 365 |
| [CppCodes/Solver_CFD_BoundaryConditions_Premixed_T.cpp](/CppCodes/Solver_CFD_BoundaryConditions_Premixed_T.cpp) | C++ | 132 | 23 | 25 | 180 |
| [CppCodes/Solver_CFD_Energy.cpp](/CppCodes/Solver_CFD_Energy.cpp) | C++ | 939 | 104 | 247 | 1,290 |
| [CppCodes/Solver_CFD_Memory.cpp](/CppCodes/Solver_CFD_Memory.cpp) | C++ | 144 | 22 | 72 | 238 |
| [CppCodes/Solver_CFD_MomentumConvection.cpp](/CppCodes/Solver_CFD_MomentumConvection.cpp) | C++ | 672 | 49 | 269 | 990 |
| [CppCodes/Solver_CFD_MomentumDiffusion.cpp](/CppCodes/Solver_CFD_MomentumDiffusion.cpp) | C++ | 776 | 81 | 115 | 972 |
| [CppCodes/Solver_CFD_PoissonCoeffs.cpp](/CppCodes/Solver_CFD_PoissonCoeffs.cpp) | C++ | 225 | 40 | 65 | 330 |
| [CppCodes/Solver_CFD_RK3_Integration.cpp](/CppCodes/Solver_CFD_RK3_Integration.cpp) | C++ | 144 | 47 | 54 | 245 |
| [CppCodes/Solver_CFD_RK4_Integration.cpp](/CppCodes/Solver_CFD_RK4_Integration.cpp) | C++ | 187 | 60 | 69 | 316 |
| [CppCodes/Solver_CFD_Utilities.cpp](/CppCodes/Solver_CFD_Utilities.cpp) | C++ | 391 | 99 | 109 | 599 |
| [CppCodes/Solver_RunSolver.cpp](/CppCodes/Solver_RunSolver.cpp) | C++ | 142 | 82 | 67 | 291 |
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

[Summary](results.md) / [Details](details.md) / [Diff Summary](diff.md) / Diff Details