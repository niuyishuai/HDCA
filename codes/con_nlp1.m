function df = con_nlp1(X, n)
x = X(1:n);
y = X(n+1:2*n);
w = X(2*n+1:3*n);
z = X(end);

% Gradient of the objective function
df = [2 * (z^2 *x - z * y) + w;
    2 * (y - z * x);
    x;
    2 * x'*(z * x - y)];

end