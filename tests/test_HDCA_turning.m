%% Test a DCA-type algorithm by turning a parameter
%% generate random data
if 1
    n=20;
    randbnd=[-1,1];
    T=rand(n)*(randbnd(2)-randbnd(1))+randbnd(1);
    mu = -min([0,eigs(T+T',1,'smallestreal')]);
    A=T+mu*eye(n);
    B=eye(n);
    x0.x = rand(n,1);
    x0.x = x0.x/sum(x0.x);
    x0.y = rand(n,1);
    x0.w = rand(n,1);
    x0.z = sum(x0.y);
end
%% Test a DCA-type algorithm by turning a parameter
modelname = 'DCP1';
algoname = 'HDCA-NI';
linesearch = 0;
nesterov = 1;
inertial = 0;
linesearch_type = 'exact';

% test solver
fprintf('Solving AEiCP via %s algorithms...\n',algoname);

% create a dc function object
dcf=dcfunc;
switch modelname
    case 'DCP1'
        dcf.f=@(X,n,A,B,opt)fobj_eval_f1(X,n,A,B,opt);
    case 'DCP2'
        dcf.f=@(X,n,A,B,opt)fobj_eval_f2(X,n,A,B,opt);
    case 'DCP3'
        dcf.f=@(X,n,A,B,opt)fobj_eval_f3(X,n,A,B,opt);
end

% create a dc problem object
mydcp=dcp(dcf,[]);

% create a dca object
mydca = dca(mydcp,x0);
mydca.A=A;
mydca.B=B;
mydca.tolf=0;
mydca.tolx=0;

mydca.verbose=0; % switch display mode
mydca.plot=false; % ploting all iterations
mydca.maxiter=200; % maxiter for DCA
mydca.savepoints=true; % save obj values of all iterations.
mydca.model=modelname;

paramlst = [0,5,10,15,20,25,30]; % parameter q for adca
%paramlst = [0.1,0.3,0.5,0.8,1,2,3,4,5]; % strong convex parameters
%paramlst = [0.01,0.03,0.05,0.08,0.1,0.5]; % strong convex parameters
%paramlst = [10,30,50,80,100,150,200,250,300,400,500]; % alphabar
%paramlst = [1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10]; % localsol_tol parameters
%paramlst = [true,false]; % inertial_conserv
valst = zeros(size(paramlst));
i=1;
for pa = paramlst
    mydca.linesearch = linesearch;
    mydca.nesterov = nesterov;
    mydca.inertial = inertial;
    mydca.strongcvx = 0.1;
    mydca.restartperiod = inf; % period for restarting nesterov
    mydca.adca_q = pa; % parameter q for adca
    mydca.localsol_tol = 1e-8;
    mydca.inertial_conserv = false;
    mydca.alphabar = 10;
    mydca.linesearch_type=linesearch_type;

    % solve model
    status=mydca.optimize();

    % get results
    hdcatime=status.time;
    hdca_times = hdcatime;
    hdca_objs = mydca.fopt;
    hdca_iters = status.iter;
    xopt=mydca.xopt(1:n);

    fprintf('Solution for %s: time %.3f sec, obj %.5e iters %d\n',algoname,hdca_times,hdca_objs,status.iter);
    % verification for AEiCP
    valst(i) = hdca_objs;
    i=i+1;
    computerr_disp(mydca.xopt,n,A,B);

    % plot figures
    figure(1);
    semilogy(status.fvallst,'LineWidth',2);
    hold all;
    figure(2);
    semilogy(status.timelst,status.fvallst,'LineWidth',2);
    hold all;
    drawnow;
end
% setup legend
cellstrcvx = cell(size(paramlst));
for i = 1:numel(paramlst)
    cellstrcvx{i} = num2str(paramlst(i));
end
figure(1)
legend(cellstrcvx,'Location','northeast');
setupfig('iter','fobj');
set(gcf, 'Position', [-2200,20,800,450]);
hold off;
figure(2)
legend(cellstrcvx,'Location','northeast');
setupfig('time (sec.)','fobj');
set(gcf, 'Position', [-1300,20,800,450]);
hold off;