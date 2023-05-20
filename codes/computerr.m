function c = computerr(x,n,A,B)
xopt = x(1:n);
lambda = (xopt'*A*xopt)/(xopt'*B*xopt);
wopt = lambda*B*xopt - A*xopt;
c = -log10(abs(xopt'*wopt) + norm(min(wopt,0)) + norm(min(xopt,0)));
end
