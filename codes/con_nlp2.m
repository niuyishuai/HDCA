function g = con_nlp2(X, n, A, B)
x = X(1:n);
y = X(n+1:2*n);
z = X(end);

% Gradient of the objective function
g = [2 * (z^2 *x - z * y) + (B+B')*x - A*y;
    2 * (y - z * x) - A'*x;
    2 * x'*(z * x - y)];

end