MagneticsSandbox
================

Magnetics simulation and theoretical analysis in matlab.

.fig
========
Fonaxis_2coilslowhysteresisGUI.fig  
Fonaxis_2coilsGUI.fig  

scripts
=========
AnalyseHessian.m	   
Analytic1D_2coils_2mags.m
computeForce.m	  
PotentialFunction_Visualize.m 
VisualizePointDipoleModel.m

functions
===========
Fonaxis_2coilsGUI.m : Works for 1D case to compute force given x1,x2 and current
Fonaxis_2coilslowhysteresisGUI.m : not working
computegradBonaxis.m   
solveforcurrent.m
checkZeroCross.m	   
computeFelectromagnet.m

--------------------------------------------------------------------------------------------------

-> /tex
1Dmathanalysis.tex : math write up for 2 electromagnets and 2 permanent magnets set up.


-> /onemagnetsim
- works for magnet dynamics in response to constant current
main_onemagnetsim.m: run this script to run simulator
compute_magnet_derivatives.m : for ode45, used in main_onemagnetsim.m
drawmagnet.m : function to represent visually a bar magnet in white|red for S|N 



