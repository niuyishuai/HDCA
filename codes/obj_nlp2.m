function f = obj_nlp2(X, n, A, B)
x = X(1:n);
y = X(n+1:2*n);
z = X(end);
% Objective function
f = norm(y - z * x)^2 + x' * (B*x - A*y);

end