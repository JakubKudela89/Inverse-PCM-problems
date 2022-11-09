%%
% This package is a MATLAB/Octave source code of CEC 2021 Bound Constrained Competition
%
%%
% AboutCEC 2021 Bound Constrained Competition:
% Ali Wagdy, Anas A Hadi, Ali K. Mohamed, Prachi Agrawal, Abhishek Kumar, and P. N. Suganthan,
% Problem Definitions and Evaluation Criteria for the CEC 2021 Special Session and Competition on Single Objective Bound Constrained Numerical Optimization,
% Technical Report, Nanyang Technological University, Singapore.
% https://www3.ntu.edu.sg/home/epnsugan/index_files/CEC2021/CEC2021-2.htm
%
%%
% Version: 0.1  Date: 2021/11/17
% Written by Anas A. Hadi (anas1401@gmail.com) and Abhishek Kumar (abhishek.kumar.eee13@iitbhu.ac.in)

clc;
clear all;

format long;
Runs=30;
NP = 100;
k=2^(-0); k=2/k;
m=1;
fhd=@(x) F2_n(x,k,m);

%optimum= [100, 1100, 700, 1900, 1700, 1600, 2100, 2200, 2400, 2500];
optimum= [0];
% C is the Parametrized Selection Vector, where:
% C1 is Shift indicator (0=0 ,1=F*)
% C2 Translation indicator (0=0, 1=oi)
% C3 Rotation indicator (0=I 1=M)

C = [1];
func=[1];

for D= [10]
    switch D
        case 10
            max_nfes=200000;
        case 20
            max_nfes=1000000;
        otherwise
            disp('Error..')
    end
    
    for i=1:size(C,1)
        fprintf('\n-------------------------------------------------------\n\n')
        run_lshade_zigzag(Runs,fhd,C(i,:),D,func,max_nfes,NP,optimum);
    end
end

