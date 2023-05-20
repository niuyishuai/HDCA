%% Test each DCA-type algorithm (load your problem data first)
% DCA|ADCA|BDCAe|BDCAa|ADCA|InDCA|HDCA-LI|HDCA-NI
% the settings for each DCA-type algorithm is:
% DCA: linesearch = 0; nesterov = 0; inertial = 0;
% BDCAe: linesearch = 1; nesterov = 0; inertial = 0; linesearch_type='exact';
% BDCAa: linesearch = 1; nesterov = 0; inertial = 0; linesearch_type='armijo';
% ADCA: linesearch = 0; nesterov = 1; inertial = 0; adca_q > 0; restartperiod = inf|>0;
% InDCA: linesearch = 0; nesterov = 0; inertial = 1;
% HDCA-LI: linesearch = 1; nesterov = 0; inertial = 1; linesearch_type='exact';
% HDCA-NI: linesearch = 0; nesterov = 1; inertial = 1; adca_q > 0; restartperiod = inf|>0;

modelname='DCP1'; % DCP1|DCP2|DCP3
algoname='HDCA-NI'; % DCA|ADCA|BDCAe|BDCAa|ADCA|InDCA|HDCA-LI|HDCA-NI
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

mydca.verbose=true; % switch display mode
mydca.plot=false; % ploting all iterations
mydca.maxiter=200; % maxiter for DCA-type methods
mydca.savepoints=true; % save obj values of all iterations.
mydca.model=modelname;

mydca.linesearch = 0;
mydca.nesterov = 1;
mydca.inertial = 1;
mydca.inertial_conserv = false;
mydca.strongcvx = 0.1;
mydca.restartperiod = inf; % period for restarting nesterov
mydca.adca_q = 10; % parameter q for adca
mydca.localsol_tol = 1e-8;
mydca.linesearch_type='exact';

% Test BDCA solver
fprintf('Solving AEiCP via HDCA algorithms...\n');

% solve model using ubdca
status=mydca.optimize();

% get results
bdcatime=status.time;
bdca_times = bdcatime;
bdca_objs = mydca.fopt;
bdca_iters = status.iter;
xopt=mydca.xopt(1:n);

fprintf('Solution for %s: time %.3f sec, obj %.5e iters %d\n',algoname,bdca_times,bdca_objs,status.iter);
% verification for EiCP
computerr_disp(mydca.xopt,n,A,B);

if mydca.savepoints
figure(1);
hold on
semilogy(status.fvallst,'LineWidth',2);
setupfig('iter','fobj');
figure(2);
hold on
semilogy(status.timelst,status.fvallst,'LineWidth',2);
setupfig('time (sec.)','fobj');
hold all;
figure(3);
hold on
plot(status.clst,'LineWidth',2);
setupfig('iter','c');
end
