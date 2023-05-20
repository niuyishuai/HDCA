function [fval,dfval] = fobj_eval_f3(X,n,A,B,opt)
x=X(1:n);
y=X(n+1:2*n);
w=X(2*n+1:3*n);
xy=x'*y;
xx=x'*x;

fval = y'*y + x'*w - xy^2/xx;
if opt~=1 % evaluate gradient
    dfval = [w - 2*(xy/xx)*y + 2*(xy^2/xx^2)*x;
        2*y - 2*(xy/xx)*x;
        x];
end
end