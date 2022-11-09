%% Main file for running the experiments
% IMPORTANT NOTE: Please note that before running this file .mex files for the simulators must be generated. 
% This can be done by running the "create_mex.m" file in the "Air-PCM_heat_exchanger_model" folder.

clear;clc;
format shortE
addpath('Air-PCM_heat_exchanger_model/');
addpath('algos_selected/PSO'); 
addpath('algos_selected/LSHADE'); 
addpath('algos_selected/DE'); 
addpath('algos_selected/CMA-ES');
addpath('algos_selected/ABC');
addpath('algos_selected/CS');
addpath('algos_selected/TLBO');


rng(1,'twister'); % seed for replicability

nruns = 10; % number of independent runs

optimum = 0; % theoretical optimum for all considered problems

prob_nr = 1; Ds = [6,6,6,6,6,7,7];  % problem dimensions for the case studies

for j=1:7       % selection of case studies to run (1-5 normal, 6-7 extended)
    D = Ds(j);
    NP = 10*D;  % population size
    max_nfes = 10000; % maximum number of function evaluations
    
switch j        % choosing the objective function corresponding to the different case studies
    case 1
        fhd = @(x) f_1(x);
    case 2
        fhd = @(x) f_2(x);
    case 3
        fhd = @(x) f_3(x);     
    case 4
        fhd = @(x) f_4(x);
    case 5
        fhd = @(x) f_5(x);
    case 6
        fhd = @(x) f_4_extended(x);
    case 7
        fhd = @(x) f_5_extended(x);
end

% running the algorithms

for i=1:nruns
fprintf('running LSHADE \n')
tic;res = run_lshade_zigzag(1,fhd,1,D,1,max_nfes,NP,optimum);display(res.bestval);time = toc;
res_LSHADE.time(i) = time;
res_LSHADE.bestval(i) = res.bestval;
res_LSHADE.progress_iters{i} = res.progress_iters;
res_LSHADE.progress_values{i} = res.progress_values;
res_LSHADE.solution(:,i) = res.solution;

fprintf('running DE \n')

tic;res = Run_DE(1,fhd,1,D,1,max_nfes,NP,optimum);display(res.bestval);time = toc;
res_DE.time(i) = time;
res_DE.bestval(i) = res.bestval;
res_DE.progress_iters{i} = res.progress_iters;
res_DE.progress_values{i} = res.progress_values;
res_DE.solution(:,i) = res.solution;
    
fprintf('running CMAES \n')

tic;res = cmaes(fhd,D,max_nfes); display(res.bestval);time = toc;
res_CMAES.time(i) = time;
res_CMAES.bestval(i) = res.bestval;
res_CMAES.progress_iters{i} = res.progress_iters;
res_CMAES.progress_values{i} = res.progress_values;
res_CMAES.solution(:,i) = res.solution;

fprintf('running PSO \n')

tic;res=PSO_func(D,NP,fhd,fix(max_nfes/NP),-100,100); display(res.bestval);time = toc;
res_PSO.time(i) = time;
res_PSO.bestval(i) = res.bestval;
res_PSO.progress_iters{i} = res.progress_iters;
res_PSO.progress_values{i} = res.progress_values;
res_PSO.solution(:,i) = res.solution;

fprintf('running ABC \n')
 
tic;res=ABC(fhd,-100*ones(1,D),100*ones(1,D),max_nfes,NP,NP*D); display(res.bestval);time = toc;
res_ABC.time(i) = time;
res_ABC.bestval(i) = res.bestval;
res_ABC.progress_iters{i} = res.progress_iters;
res_ABC.progress_values{i} = res.progress_values;
res_ABC.solution(:,i) = res.solution;

fprintf('running CS \n')
 
tic;res=cuckoo_search_new(D,fhd,-100,100,fix(max_nfes/NP),NP); display(res.bestval);time = toc;
res_CS.time(i) = time;
res_CS.bestval(i) = res.bestval;
res_CS.progress_iters{i} = res.progress_iters;
res_CS.progress_values{i} = res.progress_values;
res_CS.solution(:,i) = res.solution;


fprintf('running TLBO \n')

tic;res=SanitizedTLBO(fhd,-100*ones(1,D),100*ones(1,D),fix(max_nfes/NP),NP); display(res.bestval);time = toc;
res_TLBO.time(i) = time;
res_TLBO.bestval(i) = res.bestval;
res_TLBO.progress_iters{i} = res.progress_iters;
res_TLBO.progress_values{i} = res.progress_values;
res_TLBO.solution(:,i) = res.solution;

res_complete(prob_nr).DE = res_DE;
res_complete(prob_nr).CMAES = res_CMAES;
res_complete(prob_nr).LSHADE = res_LSHADE;
res_complete(prob_nr).PSO = res_PSO;
res_complete(prob_nr).ABC = res_ABC;
res_complete(prob_nr).CS = res_CS;
res_complete(prob_nr).TLBO = res_TLBO;
res_complete(prob_nr).D = D;
res_complete(prob_nr).f = fhd;

% saving results

switch j
    case 1
        file_name = 'results_CS1.mat';
    case 2
        file_name = 'results_CS2.mat';
    case 3
        file_name = 'results_CS3.mat';
    case 4
        file_name = 'results_CS4.mat';
    case 5
        file_name = 'results_CS5.mat';
    case 6
        file_name = 'results_CS4_extended.mat';
    case 7
        file_name = 'results_CS5_extended.mat';
end
save(file_name, 'res_complete');

end

clear res_DE res_CMAES res_LSHADE res_PSO res_ABC res_CS res_TLBO ;
prob_nr = prob_nr + 1;

end

