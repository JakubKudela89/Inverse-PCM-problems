# Inverse-PCM-problems
This repository is intended for the source code for the paper:

"Assessment of the performance of metaheuristic methods used for the inverse identification of effective heat capacity of phase change materials" (currently submitted to a journal)

It contains the MATLAB files necessary to run the experiments:

a) Simulator files for the 7 case studies (5 "normal" and 2 extended), which can be found in the folder "Air-PCM_heat_exchanger_model"

b) Codes for the selected metaheuristics, which can be found in the folder "algos_selected"

c) Main file "run_experiments.m" that can be used to replicate the experiments found in the paper

NOTE: Before running the "run_experiments.m" file, .mex files for the simulators must be created. This can be done by running the "create_mex.m" file in the "Air-PCM_heat_exchanger_model" folder.
