function [f, g] = objCon_nlp1(X, n)
x = X(1:n);
y = X(n+1:2*n);
w = X(2*n+1:3*n);
z = X(end);
% Objective function
f = norm(y - z * x)^2 + x' * w;

% Gradient of the objective function
g = [2 * (z^2 *x - z * y) + w;
    2 * (y - z * x);
    x;
    2 * x'*(z * x - y)];

end