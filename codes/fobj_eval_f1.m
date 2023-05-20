function [fval,dfval] = fobj_eval_f1(X,n,A,B,opt)
x=X(1:n);
y=X(n+1:2*n);
w=X(2*n+1:3*n);
z=X(3*n+1);
s=y-z*x;
    
if opt==1 % evaluate function only
    fval = norm(s)^2 + x'*w;
else % evaluate function and gradient
    fval = norm(s)^2 + x'*w;
    dfval = [w - 2*z*s;
        2*s;
        x;
        -2*(x'*s)];
end
end