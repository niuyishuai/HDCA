function f = obj_nlp3(X, n)
x = X(1:n);
y = X(n+1:2*n);
w = X(2*n+1:3*n);
% Objective function
f = y'*y + x' * w - (x'*y)^2/(x'*x);

end