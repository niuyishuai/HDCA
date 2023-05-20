function f= obj_nlp1(X, n)
x = X(1:n);
y = X(n+1:2*n);
w = X(2*n+1:3*n);
z = X(end);
% Objective function
f = norm(y - z * x)^2 + x' * w;
end