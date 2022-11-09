%% File for generating the mex files

clear;clc;
codegen GACR22_PRES_model_CS1 -args {[2000, 20000, 35, 23, 3, 800, 0.2]}
codegen GACR22_PRES_model_CS2 -args {[2500, 31000, 30, 15, 15, 800, 0.22]}
codegen GACR22_PRES_model_CS3 -args {[3000, 50000, 40, 5, 20, 800, 0.17]}
codegen GACR22_PRES_model_CS4 -args {[3000, 50000, 40, 5, 20, 800, 0.17]}
codegen GACR22_PRES_model_CS5 -args {[3000, 50000, 40, 5, 20, 800, 0.17]}
codegen GACR22_PRES_model_CS4_extended -args {[3000, 50000, 40, 5, 20, 800, 0.17, 2000]}
codegen GACR22_PRES_model_CS5_extended -args {[3000, 50000, 40, 5, 20, 800, 0.17, 2000]}