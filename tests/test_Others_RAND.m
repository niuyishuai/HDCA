%% Solving RAND(n) dataset using compared solvers (IPOPT, KNITRO, FILTERSD)
datasetidx=2; %RAND(datasetidx,n)
DIRNAME=sprintf('.\\RESULT_RAND%d',datasetidx);
if ~exist(DIRNAME, 'dir')
    mkdir(DIRNAME);
end
DATANAME="..\\datasets\\RAND%d\\RAND(%d,%d)_%d.mat";
RESULTNAME=".\\RESULT_RAND%d\\%s_RAND(%d,%d)_%d_%s.mat";
for n=[10,100,500]
    for modelname={'DCP1','DCP2','DCP3'}
        for nprob = 1:10
            %% Compare solvers
            filename = sprintf(DATANAME,datasetidx,datasetidx,n,nprob);
            load(filename);
            fprintf('Load data %s.........\n',filename);
            for solver = {'IPOPT','KNITRO','FILTERSD'}
                [xopt,fopt,c,iter,cputime] = othersolvers(n,A,B,x0,modelname{1},solver{1});
                save(sprintf(RESULTNAME,datasetidx,modelname{1},datasetidx,n,nprob,solver{1}),'xopt','fopt','c','iter','cputime');
            end
        end
    end
end


