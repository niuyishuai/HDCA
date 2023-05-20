%% Test HDCA-NI on NEP dataset with stopping tolerence
maxiter=1e+4;
tolf=1e-8;
datasetidx=2; %RAND(datasetidx,n)
modelname='DCP2';

algoname='HDCA-NI';
if ~exist(".\\RESULT_NEP", 'dir')
    mkdir(".\\RESULT_NEP");
end
RESULTNAME=".\\RESULT_NEP\\%s_%s_%s_tol.mat";
NEPlst=dir('..\\datasets\\NEP\\*.mat');
for ii=1:numel(NEPlst)
    filename=NEPlst(ii).name;
    load([NEPlst(ii).folder,'\',filename]);

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
    mydca.tolf=tolf;
    mydca.tolx=0;

    mydca.verbose=0; % switch display mode
    mydca.plot=false; % ploting all iterations
    mydca.maxiter=maxiter; % maxiter for DCA
    mydca.localsol_tol = 1e-8;
    mydca.savepoints=false; % save obj values of all iterations.
    mydca.model=modelname;

    mydca.linesearch = 0;
    mydca.nesterov = 1;
    mydca.inertial = 1;
    mydca.strongcvx = 0.1;
    mydca.restartperiod = inf; % period for restarting nesterov
    mydca.adca_q = 10; % parameter q for adca
    mydca.linesearch_type='exact';

    % Test HDCA-NI solver
    fprintf('Solving AEiCP via HDCA-NI algorithms...\n');

    % solve model using HDCA-NI
    status=mydca.optimize();

    % get results
    cputime=status.time;
    fopt = mydca.fopt;
    iter = status.iter;
    xopt=mydca.xopt;

    fprintf('Solution for %s using %s: time %.3f sec, obj %.5e iters %d\n',filename,algoname,cputime,fopt,iter);
    % verification for EiCP
    c = computerr_disp(xopt,n,A,B);

    save(sprintf(RESULTNAME,modelname,filename,algoname),'fopt','iter','cputime','c','xopt');

end

