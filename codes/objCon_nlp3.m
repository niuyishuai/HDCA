function [f, g] = objCon_nlp3(X, n)
x = X(1:n);
y = X(n+1:2*n);
w = X(2*n+1:3*n);
% Objective function
f = y'*y + x' * w - (x'*y)^2/(x'*x);

% Gradient of the objective function
s = (x'*y)/(x'*x);
g = [2 * s^2*x - 2*s*y + w;
    2 * y - 2*s*x;
    x];

end