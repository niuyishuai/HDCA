%% Solving NEP dataset using compared solvers (IPOPT, KNITRO, FILTERSD)
if ~exist(".\\RESULT_NEP", 'dir')
    mkdir(".\\RESULT_NEP");
end
NEPlst=dir('..\\datasets\\NEP\\*.mat');
RESULTNAME=".\\RESULT_NEP\\%s_%s_%s.mat";
for modelname={'DCP1','DCP2','DCP3'}
    for ii=1:numel(NEPlst)
        filename=NEPlst(ii).name;
        load([NEPlst(ii).folder,'\',filename]);
        fprintf('Load data %s.........\n',filename);
        for solver = {'IPOPT','KNITRO','FILTERSD'}
            [xopt,fopt,c,iter,cputime] = othersolvers(n,A,B,x0,modelname{1},solver{1});
            save(sprintf(RESULTNAME,modelname{1},filename,solver{1}),'xopt','fopt','c','iter','cputime');
        end
    end
end