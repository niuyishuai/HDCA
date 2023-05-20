%% Test HDCA-NI on RAND(n) datasets with stopping tolerence
maxiter=1e+4;
tolf=1e-8;
datasetidx=2; %RAND(datasetidx,n)
modelname='DCP1'; % model name DCP1|DCP2|DCP3

algoname='HDCA-NI';
DATANAME="..\\datasets\\RAND%d\\RAND(%d,%d)_%d.mat";
DIRNAME=sprintf('.\\RESULT_RAND%d',datasetidx);
if ~exist(DIRNAME, 'dir')
    mkdir(DIRNAME);
end
RESULTNAME=".\\RESULT_RAND%d\\%s_RAND(%d,%d)_%d_%s_tol.mat";
for n=[10,100,500]
    for nprob = 1:10
        filename = sprintf(DATANAME,datasetidx,datasetidx,n,nprob);
        load(filename);

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
        fprintf('Solving AEiCP via HDCA-NI...\n');

        % optimize model
        status=mydca.optimize();

        % get results
        cputime=status.time;
        fopt = mydca.fopt;
        iter = status.iter;
        xopt=mydca.xopt;

        fprintf('Solution for %s using %s: time %.3f sec, obj %.5e iters %d\n',filename,algoname,cputime,fopt,iter);
        % verification for EiCP
        c = computerr_disp(xopt,n,A,B);

        save(sprintf(RESULTNAME,datasetidx,modelname,datasetidx,n,nprob,algoname),'fopt','iter','cputime','c','xopt');

    end
end

