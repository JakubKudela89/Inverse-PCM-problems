clear all
format short e;
Func_Num=1;
D=10;
if D==5
        Max_FEs = 50000;
    elseif D==10
        Max_FEs = 1000000;
    elseif D==15
        Max_FEs = 3000000;
    else
        Max_FEs = 10000000;
end

k = 2^(-2); k=2/k;m=1;
fhd = @(x) F1_n(x,k,m);
Xmin=-100;
Xmax=100;
pop_size=100;
iter_max=fix(Max_FEs/pop_size);
runs=30;
for i=1:1
    func_num=Func_Num(i);
    for j=1:runs
        j,
        [gbest,gbestval,FES]= PSO_func(D,pop_size,fhd,iter_max,Xmin,Xmax);
        xbest(j,:)=gbest;
        fbest(i,j)=gbestval;
        fbest(i,j)
    end
    f_mean(i)=mean(fbest(i,:));
end

[min(fbest) median(fbest) mean(fbest) max(fbest) std(fbest)]

