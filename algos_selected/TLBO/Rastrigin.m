function F = Rastrigin(X)

[ros, ~] = size(X);
F = zeros(ros, 1);
for k = 1: ros
    x = X(k,:);
    F(k,1) = sum((x.^2 - 10.*cos(2.*pi.*x) + 10));
end