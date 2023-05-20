function [subprob,param] = initsubprob2(n,A,B,tol)
% initalize convex subproblem for DCP1
% setting parts independent to the iteration k 
% (NB: the only part depends on k is the linear coefficients, i.e., subprob.c).
%
% INPUT:
% n: size of vector x
% A,B: matrices for AEiCP
% tol: tolerence for mosek
%
% OUTPUT:
% subprob: problem structure for mosek
% param: parameters for mosek

% define some constant matrices and vectors
veconest = ones(1,n);
matid = speye(n);
vecz = zeros(n,1);
matz = sparse(2*n+7,2*n+7);

% set linear objective
%subprob.c=zeros(3*n+7,1);
% set linear constraints a
subprob.a=sparse(6+n+2,2*n+7);
% linear part for Qi
subprob.a(1,end) = -1;
subprob.a(2,end-1) = -1;
subprob.a(3,end-6) = 2;
subprob.a(3,end-5) = -1;
subprob.a(4,end-6) = -2;
subprob.a(4,end-3) = -1;
subprob.a(5,end-2) = -1;
subprob.a(6,end-4) = -1;
% Bx - Ay
subprob.a(6+1:6+n,1:2*n)=[B, -A];
subprob.a(6+n+1,1:n)=veconest;
subprob.a(6+n+2,n+1:2*n)=veconest;
subprob.a(6+n+2,2*n+1)=-1;
% set bounds
subprob.blc=[-inf;-inf;-inf;-inf;-inf;-inf;vecz;1;0];
subprob.buc=[0;0;-1;-1;0;0;inf(n,1);1;0];
subprob.blx=zeros(2*n+7,1);
subprob.bux=[];
% set Qo lower triangular only
Qo=matz;
Qo(1:n,1:n)=B+B'+matid/2;
Qo(n+1:2*n,1:n)=-(A')/2;
Qo(n+1:2*n,n+1:2*n)=2*matid + (A'*A)/2;
v1=[1 1 0 0 0 0];
v2=[0 0 1 1 0 0];
v3=[0 0 0 0 1 1];
Qo(end-5:end,end-5:end)=(v1'*v1 + v2'*v2)/8 + v3'*v3;
[subprob.qosubi,subprob.qosubj,subprob.qoval] = find(tril(Qo));
% set Qi
subprob.qcsubi=[];
subprob.qcsubj=[];
subprob.qcval=[];
subprob.qcsubk=[];
% Q1
Q1 = matz;
Q1(1:n,1:n)=2*matid;
[qcsubi,qcsubj,qcval] = find(Q1);
subprob.qcsubi = [subprob.qcsubi;qcsubi];
subprob.qcsubj = [subprob.qcsubj;qcsubj];
subprob.qcval = [subprob.qcval;qcval];
subprob.qcsubk = [subprob.qcsubk;ones(size(qcval))];
% Q2=Q3=Q4
qcsubi=2*n+1;qcsubj=2*n+1;qcval=2;
subprob.qcsubi = [subprob.qcsubi;qcsubi;qcsubi;qcsubi];
subprob.qcsubj = [subprob.qcsubj;qcsubj;qcsubj;qcsubj];
subprob.qcval = [subprob.qcval;qcval;qcval;qcval];
subprob.qcsubk = [subprob.qcsubk;2;3;4];
% Q5 lower triangular only
Q5 = matz;
Q5(1:n,1:n)=2*matid;
Q5(n+1:2*n,1:n)=2*matid;
Q5(n+1:2*n,n+1:2*n)=2*matid;
[qcsubi,qcsubj,qcval] = find(Q5);
subprob.qcsubi = [subprob.qcsubi;qcsubi];
subprob.qcsubj = [subprob.qcsubj;qcsubj];
subprob.qcval = [subprob.qcval;qcval];
subprob.qcsubk = [subprob.qcsubk;5*ones(size(qcval))];
% Q6 lower triangular only
Q6 = matz;
Q6(1:n,1:n)=2*matid;
Q6(n+1:2*n,1:n)=-2*matid;
Q6(n+1:2*n,n+1:2*n)=2*matid;
[qcsubi,qcsubj,qcval] = find(Q6);
subprob.qcsubi = [subprob.qcsubi;qcsubi];
subprob.qcsubj = [subprob.qcsubj;qcsubj];
subprob.qcval = [subprob.qcval;qcval];
subprob.qcsubk = [subprob.qcsubk;6*ones(size(qcval))];
% set initial point
%subprob.sol.itr.xx=x0;

% set stop tol, must reduce the three tols together
%param.MSK_IPAR_OPTIMIZER='MSK_OPTIMIZER_INTPNT';
%param=[];
param.MSK_DPAR_INTPNT_QO_TOL_REL_GAP = tol; 
param.MSK_DPAR_INTPNT_QO_TOL_PFEAS = tol; 
param.MSK_DPAR_INTPNT_QO_TOL_DFEAS = tol; 
end