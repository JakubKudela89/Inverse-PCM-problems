function [val] = f_3(x)
% objective function for case study 3
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

fhd = @GACR22_PRES_model_CS3_mex; opt = 0;

[m,n] = size(x); % n je dimenze

a = [1500, 10000, 10, 0.5, 0.5, 0.1]; b = [5000, 100000, 60, 50, 50, 10];

x = x/200.*(b-a)+(a+b)/2;

x(:,7) = x(:,6); x(:,6) = 800;

if m > 1 && n> 1 
    parfor i=1:m
        val(i,1) = feval(fhd,x(i,:)) - opt;
        if isnan(val(i,1))
            val(i,1) = 1e6;
        end
    end    
else
    val = feval(fhd,x) - opt;
    if isnan(val)
        val = 1e6;
    end
end

end

