%% Test all DCA-type algorithm for solving one test problem via one test DCP model
%% random dataset
if 0
    n=20;
    randbnd=[-1,1];
    T=rand(n)*(randbnd(2)-randbnd(1))+randbnd(1);
    mu = -min([0,eigs(T+T',1,'smallestreal')]);
    A=T+mu*eye(n);
    B=eye(n);
    x0.x = rand(n,1);
    x0.x = x0.x/sum(x0.x);
    x0.y = rand(n,1);
    x0.w = B*x0.x - A*x0.y;
    %x0.w = rand(n,1);
    x0.z = sum(x0.y);
    save(sprintf("RAND\\RAND(1,%d)_%s.mat",n,datestr(datetime('now'),30)),"n","A","B","x0");
end
%% Test DCA-type algorithms
% create a dc function object
dcf=dcfunc;
dcf.f=@(X,n,A,B,opt)fobj_eval_f1(X,n,A,B,opt);
%dcf.f=@(X,n,A,B,opt)fobj_eval_f2(X,n,A,B,opt);
%dcf.f=@(X,n,A,B,opt)fobj_eval_f3(X,n,A,B,opt);

% create a dc problem object
mydcp=dcp(dcf,[]);

% create a dca object
mydca = dca(mydcp,x0);
mydca.A=A;
mydca.B=B;
mydca.tolf=0;
mydca.tolx=0;

mydca.alphabar=10;
mydca.verbose=false; % switch display mode
mydca.plot=false; % ploting all iterations
mydca.maxiter=200; % maxiter for DCA
mydca.localsol_tol = 1e-8;
mydca.savepoints=true; % save obj values of all iterations.
mydca.model='DCP1'; % DCP1|DCP2|DCP3
if ~exist(".\\RESULT", 'dir')
    mkdir(".\\RESULT");
end
% Test BDCA solver
fprintf('Solving AEiCP via DCA-type algorithms...\n');
% play with linesearch, nesterov and inertial (possible to combine all of them)
% (linesearch, nesterov, inertial, 'exact'|'armijo')
testmodes={{1,0,0,'exact'},{1,0,0,'armijo'},{0,0,0,'exact'},{0,1,0,'exact'},{0,0,1,'exact'},{1,0,1,'exact'},{0,1,1,'exact'},{1,1,0,'exact'}};
%testmodes={{1,0,0,'armijo'},{0,0,0,'armijo'},{0,1,0,'armijo'},{0,0,1,'armijo'},{1,0,1,'armijo'},{0,1,1,'armijo'},{1,1,0,'armijo'}};
algoname={'BDCAe','BDCAa','DCA','ADCA','InDCA','HDCA-LI','HDCA-NI','NBDCAe'};
%algoname={'BDCAa','DCA','ADCA','InDCA','HDCA-LI','HDCA-NI','NBDCAe'};
idx=1;
for mode = testmodes
    mydca.linesearch = mode{1}{1};
    mydca.nesterov = mode{1}{2};
    mydca.inertial = mode{1}{3};
    mydca.strongcvx = 0.1;
    mydca.restartperiod = inf; % period for restarting nesterov
    mydca.adca_q = 10; % parameter q for adca
    mydca.linesearch_type=mode{1}{4};

    % solve model using ubdca
    status=mydca.optimize();

    % get results
    bdcatime=status.time;
    bdca_times = bdcatime;
    bdca_objs = mydca.fopt;
    bdca_iters = status.iter;
    xopt=mydca.xopt(1:n);

    fprintf('Solution for %s: time %.3f sec, obj %.5e iters %d\n',algoname{idx},bdca_times,bdca_objs,status.iter);
    save(sprintf("RESULT\\%s_%d_%d_%d_%s.mat",algoname{idx},mydca.linesearch,mydca.nesterov,mydca.inertial,mydca.linesearch_type),...
        "status");
    % verification for EiCP
    computerr_disp(mydca.xopt,n,A,B);
    idx = idx + 1;
end

%% Draw pictures
% compare objectives vs iterations among DC algorithms
figure(1);
clf;
for i=1:numel(algoname)
    linesearch = testmodes{i}{1};
    nesterov = testmodes{i}{2};
    inertial=testmodes{i}{3};
    linesearch_type=testmodes{i}{4};
    nn = algoname{i};
    switch nn
        case 'BDCAe'
            plotlinetype='r-';
        case 'BDCAa'
            plotlinetype='r-.';
        case 'ADCA'
            plotlinetype='b:';
        case 'DCA'
            plotlinetype='k-.';
        case 'InDCA'
            plotlinetype='b-';
        case 'HDCA-LI'
            plotlinetype='m-';
        case 'HDCA-NI'
            plotlinetype='k-';
        case 'NBDCAe'
            plotlinetype='g-';
    end

    fname = sprintf("RESULT\\%s_%d_%d_%d_%s.mat",nn,linesearch,nesterov,inertial,linesearch_type);
    load(fname);
    semilogy(status.fvallst,plotlinetype,'LineWidth',2);
    %plot(log(status.fvallst),plotlinetype,'LineWidth',2);
    hold on;
end
hold off;
legend(algoname,'Location','northeast');
setupfig('iter','fobj');
set(gcf, 'Position', [0,20,800,450]);
savefig('compDC_DCP1_fobj');

%%
% compare objectives vs time among DC algorithms
figure(2);
clf;
for i=1:numel(algoname)
    linesearch = testmodes{i}{1};
    nesterov = testmodes{i}{2};
    inertial=testmodes{i}{3};
    linesearch_type=testmodes{i}{4};
    nn = algoname{i};
    switch nn
        case 'BDCAe'
            plotlinetype='r-';
        case 'BDCAa'
            plotlinetype='r-.';
        case 'ADCA'
            plotlinetype='b:';
        case 'DCA'
            plotlinetype='k-.';
        case 'InDCA'
            plotlinetype='b-';
        case 'HDCA-LI'
            plotlinetype='m-';
        case 'HDCA-NI'
            plotlinetype='k-';
        case 'NBDCAe'
            plotlinetype='g-';
    end

    fname = sprintf("RESULT\\%s_%d_%d_%d_%s.mat",nn,linesearch,nesterov,inertial,linesearch_type);
    load(fname);
    semilogy(status.timelst,status.fvallst,plotlinetype,'LineWidth',2);
    hold on;
end
hold off;
legend(algoname,'Location','northeast');
setupfig('time (sec.)','fobj');
set(gcf, 'Position', [980,20,800,450]);
savefig('compDC_DCP1_time');
