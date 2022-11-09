function [res, X, FVAL, BestFVALIter, pop] = SanitizedTLBO(FITNESSFCN,lb,ub,T,NPop)
% Teaching Learning Based optimization (TLBO)
% SanitizedTLBO attempts to solve problems of the following forms:
%         min F(X)  subject to: lb <= X <= ub
%          X        
%                           
%  [X,FVAL,BestFVALIter, pop] = SanitizedTLBO(FITNESSFCN,lb,ub,T,NPop)
%  FITNESSFCN   - function handle of the fitness function
%  lb           - lower bounds on X
%  ub           - upper bounds on X
%  T            - number of iterations
%  NPop         - size of the population (class size)  
%  X            - minimum of the fitness function determined by SanitizedTLBO
%  FVAL         - value of the fitness function at the minima (X)
%  BestFVALIter - the best fintess function value in each iteration
%  pop          - the population at the end of the specified number of iterations 

% preallocation to store the best objective function of every iteration
% and the objective function value of every student
BestFVALIter = NaN(T,1);
obj = NaN(NPop,1);
T = T/2; %!!!two function evals per iteration
% Determining the size of the problem
D = length(lb);

% Generation of initial population
pop = repmat(lb, NPop, 1) + repmat((ub-lb),NPop,1).*rand(NPop,D);


%  Evaluation of objective function
%  Can be vectorized 
parfor p = 1:NPop
    obj(p) = FITNESSFCN(pop(p,:));
end

nr_points_to_save = 10;
RecordFEsFactor = round(logspace(log(2*NPop)/log(10),log(T*NPop)/log(10),nr_points_to_save));
Run_RecordFEsFactor=RecordFEsFactor;
current_eval = NPop; run_funcvals = []; run_funccalls = [];

run_funcvals = [run_funcvals; min(obj)];
        run_funccalls = [run_funccalls; current_eval];
        fprintf('fes: %u, val: %5.3e\n', current_eval, min(obj))


for gen = 1: T
    
    % Partner selection for all students
    % Note that randperm has been used to speedup the partner selection.
    Partner = randperm(NPop);
    % There is a remote possibility that the ith student will have itself as its partner
    % No experiment is available in literature on the disadvantages of
    % a solution having itself as partner solution.
    pop_old = pop; obj_old = obj;
    parfor i = 1:NPop
        
        % ----------------Begining of the Teacher Phase for ith student-------------- %
        mean_stud = mean(pop_old);
        
        % Determination of teacher
        [~,ind] = min(obj_old);
        best_stud = pop_old(ind,:);
        
        % Determination of the teaching factor
        TF = randi([1 2],1,1);
        
        % Generation of a new solution
        NewSol = pop_old(i,:) + rand(1,D).*(best_stud - TF*mean_stud);
        
        % Bounding of the solution
        NewSol = max(min(ub, NewSol),lb);
        
        % Evaluation of objective function
        NewSolObj = FITNESSFCN(NewSol);
        
        % Greedy selection
        if (NewSolObj < obj_old(i))
            pop(i,:) = NewSol;
            obj(i) = NewSolObj;
        end
        % ----------------Ending of the Teacher Phase for ith student-------------- %
        
        
        % ----------------Begining of the Learner Phase for ith student-------------- %
        % Generation of a new solution
        if (obj(i)< obj_old(Partner(i)))
            NewSol = pop_old(i,:) + rand(1, D).*(pop_old(i,:)- pop_old(Partner(i),:));
        else
            NewSol = pop_old(i,:) + rand(1, D).*(pop_old(Partner(i),:)- pop_old(i,:));
        end
        
        % Bounding of the solution
        NewSol = max(min(ub, NewSol),lb);
        
        % Evaluation of objective function
        NewSolObj =  FITNESSFCN(NewSol);
        
        % Greedy selection
        if(NewSolObj< obj(i))
            pop(i,:) = NewSol;
            obj(i) = NewSolObj;
        end
        % ----------------Ending of the Learner Phase for ith student-------------- %
        
    end
    
    % This is not part of the algorithm but is used to keep track of the
    % best solution determined till the current iteration
    [BestFVALIter(gen),ind] = min(obj);
    current_eval = current_eval + 2*NPop;
    
    run_funcvals = [run_funcvals; BestFVALIter(gen)];
        run_funccalls = [run_funccalls; current_eval];
        fprintf('fes: %u, val: %5.3e\n', current_eval, BestFVALIter(gen))
end

% Extracting the best solution
X = pop(ind,:);
FVAL = BestFVALIter(gen);

run_funcvals(end) = FVAL;
res.solution = X;
res.bestval = FVAL;
res.progress_iters = run_funccalls;
res.progress_values = run_funcvals;