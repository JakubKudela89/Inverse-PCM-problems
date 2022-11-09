clc
clear
close all


rng(5,'twister')

FITNESSFCN = @Rastrigin;

D = 2;
lb = -5.12*ones(1,D);
ub = 5.12*ones(1,D);

NumberFoods = 10; % Population Size
MaxFe = 4000; % Maximum number of functional evaluations
limit = NumberFoods*D; % NumberFoods*D

[X, FVAL, Foods, ObjVal, BestFVALCycle] = ABC(FITNESSFCN,lb,ub,MaxFe,NumberFoods,limit);

display(['The minimum point is ', num2str(X)])
display(['The fitness function value at the mimimum point is ', num2str(FVAL)])
display(['The number of cycles utilized is ', num2str(length(BestFVALCycle))])

plot(BestFVALCycle,'r*')
xlabel('Cycle Number')
ylabel('Value of Fitness function')
grid on