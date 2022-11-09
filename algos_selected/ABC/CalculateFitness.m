function Fitness = CalculateFitness(ObjVal)

% Note that minimizing the objective function is equivalent to maximization
% of fitness as they are inversely related.

% Fitness = zeros(size(ObjVal));
% Equation 3 of Reference [1]

% Determine the fitness of the solutions whose objective function is
% greater than or equal to 0
Fitness(ObjVal >= 0) = 1./(ObjVal(ObjVal >= 0)+1);

% Determine the fitness of the solutions whose objective function is
% less than 0
Fitness(ObjVal<0) = 1+abs(ObjVal(ObjVal<0));
