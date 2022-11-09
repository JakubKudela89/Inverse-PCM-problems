function res = cmaes(fhd,D,max_nfes)
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA108
% Project Title: Covariance Matrix Adaptation Evolution Strategy (CMA-ES)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

%% Problem Settings

CostFunction=fhd;   % Cost Function

nVar=D;                % Number of Unknown (Decision) Variables

VarSize=[1 nVar];       % Decision Variables Matrix Size

VarMin=-100;             % Lower Bound of Decision Variables
VarMax= 100;             % Upper Bound of Decision Variables

%% CMA-ES Settings

% Population Size (and Number of Offsprings)
lambda=(4+round(3*log(nVar)))*10;
MaxIt = fix(max_nfes/(lambda+1));

% Number of Parents
mu=round(lambda/2);

% Parent Weights
w=log(mu+0.5)-log(1:mu);
w=w/sum(w);

% Number of Effective Solutions
mu_eff=1/sum(w.^2);

% Step Size Control Parameters (c_sigma and d_sigma);
sigma0=0.3*(VarMax-VarMin);
cs=(mu_eff+2)/(nVar+mu_eff+5);
ds=1+cs+2*max(sqrt((mu_eff-1)/(nVar+1))-1,0);
ENN=sqrt(nVar)*(1-1/(4*nVar)+1/(21*nVar^2));

% Covariance Update Parameters
cc=(4+mu_eff/nVar)/(4+nVar+2*mu_eff/nVar);
c1=2/((nVar+1.3)^2+mu_eff);
alpha_mu=2;
cmu=min(1-c1,alpha_mu*(mu_eff-2+1/mu_eff)/((nVar+2)^2+alpha_mu*mu_eff/2));
hth=(1.4+2/(nVar+1))*ENN;

%% Initialization
run_funcvals = [];
run_funccalls = [];

ps=cell(MaxIt,1);
pc=cell(MaxIt,1);
C=cell(MaxIt,1);
sigma=cell(MaxIt,1);

ps{1}=zeros(VarSize);
pc{1}=zeros(VarSize);
C{1}=eye(nVar);
sigma{1}=sigma0;

empty_individual.Position=[];
empty_individual.Step=[];
empty_individual.Cost=[];

M=repmat(empty_individual,MaxIt,1);
M(1).Position=unifrnd(VarMin,VarMax,VarSize);
M(1).Step=zeros(VarSize);
M(1).Cost=CostFunction(M(1).Position);

BestSol=M(1);

BestCost=zeros(MaxIt,1);



%% CMA-ES Main Loop

nr_points_to_save = 10;
RecordFEsFactor = round(logspace(log(2*lambda)/log(10),log(max_nfes)/log(10),nr_points_to_save));
Run_RecordFEsFactor=RecordFEsFactor;
current_eval = 0; 
for g=1:MaxIt
    
    % Generate Samples
    pop=repmat(empty_individual,lambda,1);
    for i=1:lambda
        try
        pop(i).Step=mvnrnd(zeros(VarSize),C{g});
        catch
           C{g} = C{g} + 1e-15*eye(D);
           pop(i).Step=mvnrnd(zeros(VarSize),C{g});
        end
        pop(i).Position=M(g).Position+sigma{g}*pop(i).Step;
    end
    parfor i=1:lambda
        pop(i).Cost=CostFunction(pop(i).Position);
    end

   for i=1:lambda         
        % Update Best Solution Ever Found
        if pop(i).Cost<BestSol.Cost
            BestSol=pop(i);
        end
   end
    current_eval = current_eval + lambda;
    while(current_eval>=Run_RecordFEsFactor(1))
                    %run_funcvals = [run_funcvals; BestSol.Cost];
                    %run_iters = [run_iters; current_eval];
                    Run_RecordFEsFactor(1)=[];
                    if isempty(Run_RecordFEsFactor)
                        break
                    end
    end
    run_funcvals = [run_funcvals; BestSol.Cost];
    run_funccalls = [run_funccalls; current_eval];
    fprintf('fes: %u, val: %5.3e\n', current_eval, BestSol.Cost)
    % Sort Population
    Costs=[pop.Cost];
    [Costs, SortOrder]=sort(Costs);
    pop=pop(SortOrder);
  
    % Save Results
    BestCost(g)=BestSol.Cost;
    if BestCost(g) <= 1e-8
        BestSol.Cost = 0;
        BestCost(g)= 0;
       break 
    end
    % Display Results
    
    
    % Exit At Last Iteration
    if g==MaxIt
        break;
    end
        
    % Update Mean
    M(g+1).Step=0;
    for j=1:mu
        M(g+1).Step=M(g+1).Step+w(j)*pop(j).Step;
    end
    M(g+1).Position=M(g).Position+sigma{g}*M(g+1).Step;
    M(g+1).Cost=CostFunction(M(g+1).Position);
    if M(g+1).Cost<BestSol.Cost
        BestSol=M(g+1);
    end
    
    % Update Step Size
    try
    ps{g+1}=(1-cs)*ps{g}+sqrt(cs*(2-cs)*mu_eff)*M(g+1).Step/chol(C{g}+1e-6*eye(D))';
    catch
       v = 0; 
    end
    sigma{g+1}=sigma{g}*exp(cs/ds*(norm(ps{g+1})/ENN-1))^0.3;
    
    % Update Covariance Matrix
    if norm(ps{g+1})/sqrt(1-(1-cs)^(2*(g+1)))<hth
        hs=1;
    else
        hs=0;
    end
    delta=(1-hs)*cc*(2-cc);
    pc{g+1}=(1-cc)*pc{g}+hs*sqrt(cc*(2-cc)*mu_eff)*M(g+1).Step;
    C{g+1}=real((1-c1-cmu)*C{g}+c1*(pc{g+1}'*pc{g+1}+delta*C{g}));
    for j=1:mu
        C{g+1}=real(C{g+1}+cmu*w(j)*pop(j).Step'*pop(j).Step);
    end
    
    % If Covariance Matrix is not Positive Defenite or Near Singular
    [V, E]=eig(C{g+1});
    if any(diag(E)<0)
        E=max(E,0);
        C{g+1}=real(V*E/V);
    end
    C{g+1} = 0.5*(C{g+1}+C{g+1}');
end

%run_funcvals(end) = BestSol.Cost;
res.solution = BestSol.Position;
res.bestval = BestSol.Cost;
res.progress_iters = run_funccalls;
res.progress_values = run_funcvals;

%% Display Results
%disp(['Iteration ' num2str(g) ': Best Cost = ' num2str(BestCost(g))]);
% figure;
% % plot(BestCost, 'LineWidth', 2);
% semilogy(BestCost, 'LineWidth', 2);
% xlabel('Iteration');
% ylabel('Best Cost');
% grid on;
end
