function [fval,dfval] = fobj_eval_f2(X,n,A,B,opt)
% opt=1: compute fval only
% opt=0: compute both fval and dfval
x=X(1:n);
y=X(n+1:2*n);
z=X(2*n+1);
s=y-z*x;

if opt==1 % evaluate function only
    fval = norm(s)^2 + x'*(B*x-A*y);
else % evaluate function and gradient
    fval = norm(s)^2 + x'*(B*x-A*y);
    dfval = [(B+B')*x - 2*z*s - A*y;
        2*s - A'*x;
        -2*(x'*s)];
end
end