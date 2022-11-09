function [temp1,temp2,temp3,temp4] = GenerateNewSol(FITNESSFCN, lb, ub, NumberFoods, i, Foods, Fitness, trial, ObjVal, D)

% Randomly select the variable that is to be changed
% Reference [1]; Page 122; index j
j = randi(D,1);

% Randomly select the neighbour
% Reference [1]; Page 122; index k
k = randi(NumberFoods,1);

% Ensuring that the neighbour is different from the current solution i
while (k == i)
    k = randi(NumberFoods,1);
end

sol = Foods(i,:);

% Generating a random number between -1 and 1
Phi = -1 + (1-(-1))*rand;
% Generating a new solution using equation (2) of Reference [1]
sol(j) = Foods(i,j)+ Phi*(Foods(i,j) - Foods(k,j));

% if the variable has a value outside its boundaries, clipping it to the
% boundary; Reference [1] Page 122
sol(j) = max(sol(j),lb(j));
sol(j) = min(sol(j),ub(j));

% Determining the objective function and fitness value
ObjValSol = FITNESSFCN(sol);
FitnessSol = CalculateFitness(ObjValSol);
%NumEval = NumEval + 1;


% Note that minimizing the objective function is equivalent to maximization
% of fitness as they are inversely related.

% Apply a greedy selection between the current solution i and the newly generated solution
if (FitnessSol > Fitness(i)) 
    % The new solution enters the pool of solutions as it has a better
    % FITNESS function value. 
    Foods(i,:) = sol;
    Fitness(i) = FitnessSol;
    ObjVal(i) = ObjValSol;
    trial(i) = 0;
else
    % Increase the value of the trial counter
    trial(i) = trial(i)+1; 
end
temp1 = trial(i);temp2 = Foods(i,:);temp3 = Fitness(i);temp4 = ObjVal(i);