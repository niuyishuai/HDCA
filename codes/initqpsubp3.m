function P = initqpsubp3(n,eta,rho,A,B)
matid = eye(n);
matz = zeros(n);

P.Q = rho*eye(3*n) + [(eta+1/2)*matid, matz, matid/2; matz, (2+eta)*matid, matz; matid/2, matz, matid/2];
% P.c need to be computed later
P.Aeq =[B, -A, -matid; ones(1,n), zeros(1,2*n)]; 
P.beq = [zeros(n,1);-1];
P.Aineq = eye(3*n);
P.bineq = zeros(3*n,1);
end