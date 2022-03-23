1. The model is presented in this file. model2D.m is the main function which is the entrance of the program that call the supporting function in sequence. PsiDef.m defines the stream function (Psi). VelCalc.m calculates the SUDM velocity field. WeightCalc.m combines the coefficients from VelValc.m and initCon initiate the concentration and the fractionation coefficients. Finally, AdvDiffGPU.m finish the iteration and reach a steady state solution.

2. The isopycnal surface with photosynthesis model σ<sub>θ</sub> 25.8-26.2 is in the folder of 25.8-26.2. There are instructions in model2D on how to switch to respiration-only mode.

3. The isopycnal surface σ<sub>θ</sub> 26.5-26.9 is in the folder of 26.5-26.9. There are instructions in model2D on how to switch to implicit-photosynthesis mode.
