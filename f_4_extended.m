function [val] = f_4_extended(x)
% objective function for case study 4 with the extended model
if nargin == 0
    val.nx = 0;
    val.ng = 0;
    val.nh = 0;
    val.xl = @(i) -100;
    val.xu = @(i) 100;
    val.fmin = @(i) 0;
    val.xmin = @(i) 0;
    return
end

fhd = @GACR22_PRES_model_CS4_extended_mex; opt = 0;

[m,n] = size(x); % n je dimenze

a = [1500, 10000, 10, 0.5, 0.5, 0.1,1500]; b = [5000, 100000, 60, 50, 50, 10, 5000];

x = x/200.*(b-a)+(a+b)/2;

x(:,7:8) = x(:,6:7); x(:,6) = 800;

if m > 1 && n> 1 
    parfor i=1:m
        val(i,1) = feval(fhd,x(i,:)) - opt;
        if isnan(val(i,1))
            val(i,1) = 1e5;
        end
    end    
else
    val = feval(fhd,x) - opt;
    if isnan(val)
        val = 1e5;
    end
end

end

