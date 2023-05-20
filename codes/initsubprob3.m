function [subprob,param] = initsubprob3(n,A,B,tol)
% initalize convex subproblem for DCP3
% setting parts independent to the iteration k 
% (NB: the only part depends on k is the linear coefficients, i.e., subprob.c).
%
% INPUT:
% n: size of vector x
% A,B: matrices for AEiCP
% eta: parameter for DC decomposition
% tol: tolerence for mosek
%
% OUTPUT:
% subprob: problem structure for mosek
% param: parameters for mosek

% define some constant matrices and vectors
veconest = ones(1,n);
matid = speye(n);
vecz = zeros(n,1);


% set linear objective
subprob.c=[vecz;-ones(n,1);vecz];

% set linear constraints a
subprob.a=sparse([B, -A, -matid; veconest, vecz', vecz']);
%subprob.a=sparse(n+1,3*n);
%subprob.a(1:n,1:3*n)=[B, -A, -matid];
%subprob.a(n+1,1:n)=veconest;
% set bounds
subprob.blc=[vecz;1];
subprob.buc=[vecz;1];
subprob.blx=zeros(3*n,1);
subprob.bux=[];
% set initial point
%subprob.sol.itr.xx=x0;


% set stop tol, must reduce the three tols together
param.MSK_DPAR_INTPNT_QO_TOL_REL_GAP = tol; 
param.MSK_DPAR_INTPNT_QO_TOL_PFEAS = tol; 
param.MSK_DPAR_INTPNT_QO_TOL_DFEAS = tol; 
end