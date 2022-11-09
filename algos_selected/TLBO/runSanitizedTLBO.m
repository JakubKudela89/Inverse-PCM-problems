clc
clear 
close all


rng(2,'twister')

FITNESSFCN = @Rastrigin;

lb = -5.12*ones(1,2);
ub = 5.12*ones(1,2);

NPop = 50;
T = 90;

[X,FVAL,BestFVALIter] = SanitizedTLBO(FITNESSFCN,lb,ub,T,NPop);

 display(['The minimum point is ', num2str(X)])
 display(['The fitness function value at the mimimum point is ', num2str(FVAL)])
 D = length(lb);
 display(['The number of fitness function evaluation is ', num2str(NPop+2*NPop*T)])
 
 plot(1:T,BestFVALIter,'r*')
 xlabel('Iteration Number')
 ylabel('Value of Fitness function')
 grid on
 
     
