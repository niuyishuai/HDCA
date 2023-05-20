NEPlst=dir('NEP\\*.mtx');
for i=1:length(NEPlst)
    fname=NEPlst(i).name;
    fprintf('Loading model from %s\n',fname);
    [Aorg,n]=mmread(sprintf('NEP\\%s',fname));
    AA=(Aorg+Aorg');
    B=eye(n);
    mu = abs(min([eigs(AA, 1, 'smallestreal'),0]))+1;
    % compute best mu by solving an SDP
    %varmu=sdpvar(1);
    %optimize(A+varmu*B>=0,varmu,sdpsettings('solver','mosek','verbose',0));
    %mu=value(varmu)+1;
    A = Aorg + mu*B;
    condA = cond(Aorg);
    try 
        chol(A);
    catch
        error('non PD matrix A');
    end
    x0.x = rand(n,1);
    x0.x = x0.x/sum(x0.x);
    x0.y = rand(n,1);
    x0.w = B*x0.x - A*x0.y;
    x0.z = sum(x0.y);
    save(sprintf('NEP\\NEP_%s.mat',fname),'n','A','B','Aorg','condA','mu','x0');
end