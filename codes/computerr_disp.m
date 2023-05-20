function c = computerr_disp(x,n,A,B)
xopt = x(1:n);
lambda = (xopt'*A*xopt)/(xopt'*B*xopt);
wopt = lambda*B*xopt - A*xopt;
%zopt = 1/lambda;
%yopt = zopt*xopt;
c = -log10(abs(xopt'*wopt) + norm(min(wopt,0)) + norm(min(xopt,0)));
%norm(min(wopt,0))
fprintf(' * <x, w> :             %e\n',xopt'*wopt);
fprintf(' * w positivity error : %e\n',norm(min(wopt,0)));
fprintf(' * x positivity error : %e\n',norm(min(xopt,0)));
fprintf(' * lambda :             %f\n',lambda);
fprintf(' * c :                  %d\n',c);
end