%% Test all DCA-type algorithms for solving RAND(n) datasets
maxiter=200;
datasetidx=2; %RAND(datasetidx,n)
modelname='DCP1';
% play with linesearch, nesterov and inertial (possible to combine all of them)
% (linesearch, nesterov, inertial, 'exact'|'armijo')
testmodes={{0,0,0,'exact'},{1,0,0,'exact'},{1,0,0,'armijo'},{0,1,0,'exact'},{0,0,1,'exact'},{1,0,1,'exact'},{0,1,1,'exact'}};
%testmodes={{0,0,0,'exact'},{1,0,0,'exact'},{1,0,0,'armijo'},{0,1,0,'exact'},{0,0,1,'exact'},{1,0,1,'exact'},{0,1,1,'exact'},{1,1,0,'exact'}};
algoname={'DCA','BDCAe','BDCAa','ADCA','InDCA','HDCA-LI','HDCA-NI'};
DATANAME="..\\datasets\\RAND%d\\RAND(%d,%d)_%d.mat";
DIRNAME=sprintf('.\\RESULT_RAND%d',datasetidx);
if ~exist(DIRNAME, 'dir')
    mkdir(DIRNAME);
end
RESULTNAME=".\\RESULT_RAND%d\\%s_RAND(%d,%d)_%d_%s.mat";
PLOTRESULTNAME=".\\RESULT_RAND%d\\%s_PLOTRAND(%d,%d).mat";
for n=[10,100,500]
    fvallst = zeros(numel(testmodes),maxiter);
    timelst = fvallst;
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

        % Test BDCA solver
        fprintf('Solving AEiCP via DCA-type algorithms...\n');

        % create a dca object
        mydca = dca(mydcp,x0);
        mydca.A=A;
        mydca.B=B;
        mydca.tolf=0;
        mydca.tolx=0;

        mydca.verbose=false; % switch display mode
        mydca.plot=false; % ploting all iterations
        mydca.maxiter=maxiter; % maxiter for DCA
        mydca.localsol_tol = 1e-8;
        mydca.savepoints=true; % save obj values of all iterations.
        mydca.model=modelname;


        idx = 1; % algorithm index
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

            save(sprintf(RESULTNAME,datasetidx,modelname,datasetidx,n,nprob,algoname{idx}),'status','mydca','maxiter');

            fprintf('Solution for %s using %s: time %.3f sec, obj %.5e iters %d\n',filename,algoname{idx},bdca_times,bdca_objs,status.iter);
            % verification for EiCP
            computerr_disp(mydca.xopt,n,A,B);

            fvallst(idx,:) = fvallst(idx,:) + status.fvallst';
            timelst(idx,:) = timelst(idx,:) + status.timelst';
            idx = idx + 1;
        end
    end

    fvallst = fvallst/10;
    timelst = timelst/10;
    save(sprintf(PLOTRESULTNAME,datasetidx,modelname,datasetidx,n),"timelst","fvallst");

    %% Draw pictures
    % compare objectives vs iterations among DC algorithms
    figure(1);
    clf;
    figure(2);
    clf;
    for i=1:numel(algoname)
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
        figure(1);
        semilogy(fvallst(i,:),plotlinetype,'LineWidth',2);
        hold all;
        figure(2);
        semilogy(timelst(i,:),fvallst(i,:),plotlinetype,'LineWidth',2);
        hold all;
    end
    figure(1);
    legend(algoname,'Location','best');
    setupfig('iter','Average fobj');
    set(gcf, 'Position', [-1800,20,800,450]);
    hold off;
    figure(2);
    legend(algoname,'Location','northeast');
    setupfig('cpu time (sec.)','Average fobj');
    set(gcf, 'Position', [-980,20,800,450]);
    hold off;
    % bar plot
    figure(3);
    bardata = timelst(1:numel(algoname),end);
    bar(bardata);
    xticklabels(algoname);

    axes1 = gca;
    set(axes1,'looseInset',[0 0 0 0]);
    fontsz = 13;

    set(axes1,'FontWeight','bold','FontSize',fontsz,'LineWidth',1);

    ylabel('Average cpu time (sec.)','FontWeight','bold','FontSize',fontsz);

    set(gcf,'InvertHardcopy','off','PaperUnits','points',...
        'Color',[1 1 1],...
        'Renderer','painters',...
        'position',[100 300 400 450]);
    drawnow;
end
