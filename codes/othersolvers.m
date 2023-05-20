function [xopt,fopt,c,iters,cputime] = othersolvers(n,A,B,x0,modelname,solvername)
% [xopt,fopt,iters,cputime] = othersolvers(n,A,B,x0,modelname,solvername)
% 
% modelname: DCP1|DCP2|DCP3
% solvername: KNITRO|FILTERSD|FMINCON|IPOPT
modelname = upper(modelname);
solvername = upper(solvername);

switch solvername
    case 'KNITRO'
        fprintf('Solving model %s using KNITRO\n',modelname);
        switch modelname
            case 'DCP1'
                x0knitro = [x0.x;x0.y;x0.w;x0.z];
                fobj = @(x) objCon_nlp1(x,n);
                % Set variable bounds
                lb = zeros(3 * n+1, 1);
                ub = inf(3 * n + 1, 1);

                % Set linear constraints
                Aeq = [B, -A, -eye(n), zeros(n, 1);
                    ones(1, n), zeros(1, n), zeros(1, n), 0;
                    zeros(1, n), ones(1, n), zeros(1, n), -1];
                beq = [zeros(n, 1); 1; 0];

                Aineq = [];
                bineq = [];
            case 'DCP2'
                x0knitro = [x0.x;x0.y;x0.z];
                fobj = @(x) objCon_nlp2(x,n,A,B);
                % Set variable bounds
                lb = zeros(2 * n+1, 1);
                ub = inf(2 * n + 1, 1);

                % Set linear constraints
                Aeq = [ones(1, n), zeros(1, n), 0;
                    zeros(1, n), ones(1, n), -1];
                beq = [1; 0];

                Aineq = [-B, A, zeros(n, 1)];
                bineq = zeros(n,1);
            case 'DCP3'
                x0knitro = [x0.x;x0.y;x0.w];
                fobj = @(x) objCon_nlp3(x,n);
                % Set variable bounds
                lb = zeros(3 * n, 1);
                ub = inf(3 * n, 1);

                % Set linear constraints
                Aeq = [B, -A, -eye(n);
                    ones(1, n), zeros(1, n), zeros(1, n)];
                beq = [zeros(n, 1); 1];

                Aineq = [];
                bineq = [];
        end
        option = optimset('Display','off','GradObj','on');
        % 'MaxIter',1000,'TolX', 1e-15, 'TolFun', 1e-8, 'TolCon', 1e-8);

        tic
        [xopt,fopt,exitflag,info]=knitromatlab(fobj,x0knitro,Aineq,bineq,Aeq,beq,lb,ub,[],[],option);

        % get results
        cputime=toc;
        iters = info.iterations;
        fprintf('Solution of KNITRO for model %s: time %.3f sec, obj %.5e iters %d\n',modelname,cputime,fopt,iters);
        c=computerr_disp(xopt,n,A,B);
        %% Test FILTERSD
    case 'FILTERSD'
        fprintf('Solving model %s using FLISTERSD\n',modelname);
        switch modelname
            case 'DCP1'
                x0filtersd = [x0.x;x0.y;x0.w;x0.z];
                fobj = @(x) obj_nlp1(x, n); % output first term only
                gradfobj = @(x) con_nlp1(x, n); % output the second term only
                % Set variable bounds
                lb = zeros(3 * n+1, 1);
                ub = inf(3 * n + 1, 1);

                % Set linear constraints
                Aeq = [B, -A, -eye(n), zeros(n, 1);
                    ones(1, n), zeros(1, n), zeros(1, n), 0;
                    zeros(1, n), ones(1, n), zeros(1, n), -1];
                beq = [zeros(n, 1); 1; 0];
                lcon=@(x) Aeq*x - beq;
                ljac=@(x) Aeq;
                cl=zeros(n+2,1);
                cu=zeros(n+2,1);
            case 'DCP2'
                x0filtersd = [x0.x;x0.y;x0.z];
                fobj = @(x) obj_nlp2(x, n, A, B); % output first term only
                gradfobj = @(x) con_nlp2(x, n, A, B); % output the second term only
                % Set variable bounds
                lb = zeros(2 * n+1, 1);
                ub = inf(2 * n + 1, 1);

                % Set linear constraints
                Aeq = [ones(1, n), zeros(1, n), 0;
                    zeros(1, n), ones(1, n), -1];
                beq = [1; 0];

                Aineq = [-B, A, zeros(n, 1)];
                bineq = zeros(n,1);
                AA = [Aeq;Aineq];
                bb = [beq;bineq];
                lcon=@(x) AA*x - bb;
                ljac=@(x) AA;
                cl=[0;0;-inf(n,1)];
                cu=[0;0;zeros(n,1)];
            case 'DCP3'
                x0filtersd = [x0.x;x0.y;x0.w];
                fobj = @(x) obj_nlp3(x, n); % output first term only
                gradfobj = @(x) con_nlp3(x, n); % output the second term only
                % Set variable bounds
                lb = zeros(3 * n, 1);
                ub = inf(3 * n, 1);

                % Set linear constraints
                Aeq = [B, -A, -eye(n);
                    ones(1, n), zeros(1, n), zeros(1, n)];
                beq = [zeros(n, 1); 1];
                lcon=@(x) Aeq*x - beq;
                ljac=@(x) Aeq;
                cl=zeros(n+1,1);
                cu=zeros(n+1,1);
        end

        tic
        [xopt,fopt,exitflag,info] = filtersd(fobj,gradfobj,x0filtersd,lb,ub,lcon,ljac,cl,cu);

        % get results
        cputime=toc;
        iters = info.niter;
        fprintf('Solution of FILTERSD for model %s: time %.3f sec, obj %.5e iters %d\n',modelname,cputime,fopt,iters);
        c=computerr_disp(xopt,n,A,B);

        %% Test FMINCON
    case 'FMINCON'
        fprintf('Solving model %s using FMINCON\n',modelname);
        switch modelname
            case 'DCP1'
                x0fmincon = [x0.x;x0.y;x0.w;x0.z];
                fobj = @(x) objCon_nlp1(x,n);
                % Set variable bounds
                lb = zeros(3 * n+1, 1);
                ub = inf(3 * n + 1, 1);

                % Set linear constraints
                Aeq = [B, -A, -eye(n), zeros(n, 1);
                    ones(1, n), zeros(1, n), zeros(1, n), 0;
                    zeros(1, n), ones(1, n), zeros(1, n), -1];
                beq = [zeros(n, 1); 1; 0];

                Aineq = [];
                bineq = [];
            case 'DCP2'
                x0fmincon = [x0.x;x0.y;x0.z];
                fobj = @(x) objCon_nlp2(x,n,A,B);
                % Set variable bounds
                lb = zeros(2 * n+1, 1);
                ub = inf(2 * n + 1, 1);

                % Set linear constraints
                Aeq = [ones(1, n), zeros(1, n), 0;
                    zeros(1, n), ones(1, n), -1];
                beq = [1; 0];

                Aineq = [-B, A, zeros(n, 1)];
                bineq = zeros(n,1);
            case 'DCP3'
                x0fmincon = [x0.x;x0.y;x0.w];
                fobj = @(x) objCon_nlp3(x,n);
                % Set variable bounds
                lb = zeros(3 * n, 1);
                ub = inf(3 * n, 1);

                % Set linear constraints
                Aeq = [B, -A, -eye(n);
                    ones(1, n), zeros(1, n), zeros(1, n)];
                beq = [zeros(n, 1); 1];

                Aineq = [];
                bineq = [];
        end
        option = optimset('Display','off','GradObj','on');
        % 'MaxIter',1000,'TolX', 1e-15, 'TolFun', 1e-8, 'TolCon', 1e-8);

        tic
        [xopt,fopt,exitflag,info]=fmincon(fobj,x0fmincon,Aineq,bineq,Aeq,beq,lb,ub,[],option);

        % get results
        cputime=toc;
        iters = info.iterations;
        fprintf('Solution of FMINCON for model %s: time %.3f sec, obj %.5e iters %d\n',modelname,cputime,fopt,iters);
        c=computerr_disp(xopt,n,A,B);

        %% Test IPOPT
    case  'IPOPT'
        fprintf('Solving model %s using IPOPT\n',modelname);
        switch modelname
            case 'DCP1'
                x0ipopt = [x0.x;x0.y;x0.w;x0.z];
                funchandle.objective = @(x) obj_nlp1(x, n); % output first term only
                funchandle.gradient = @(x) con_nlp1(x, n); % output the second term only
                % Set variable bounds
                option.lb = zeros(3 * n+1, 1);
                option.ub = inf(3 * n + 1, 1);

                % Set linear constraints
                Aeq = [B, -A, -eye(n), zeros(n, 1);
                    ones(1, n), zeros(1, n), zeros(1, n), 0;
                    zeros(1, n), ones(1, n), zeros(1, n), -1];
                beq = [zeros(n, 1); 1; 0];
                funchandle.constraints=@(x) Aeq*x - beq;
                funchandle.jacobian=@(x) sparse(Aeq);
                [Aeq_i, Aeq_j] = find(Aeq);
                funchandle.jacobianstructure = @() sparse(Aeq_i, Aeq_j, ones(length(Aeq_i), 1), size(Aeq, 1), size(Aeq, 2));
                option.cl=zeros(n+2,1);
                option.cu=zeros(n+2,1);
            case 'DCP2'
                x0ipopt = [x0.x;x0.y;x0.z];
                funchandle.objective = @(x) obj_nlp2(x, n, A, B); % output first term only
                funchandle.gradient = @(x) con_nlp2(x, n, A, B); % output the second term only
                % Set variable bounds
                option.lb = zeros(2 * n+1, 1);
                option.ub = inf(2 * n + 1, 1);

                % Set linear constraints
                Aeq = [ones(1, n), zeros(1, n), 0;
                    zeros(1, n), ones(1, n), -1];
                beq = [1; 0];

                Aineq = [-B, A, zeros(n, 1)];
                bineq = zeros(n,1);
                AA = [Aeq;Aineq];
                bb = [beq;bineq];
                funchandle.constraints=@(x) AA*x - bb;
                funchandle.jacobian=@(x) sparse(AA);
                [AA_i, AA_j] = find(AA);
                funchandle.jacobianstructure = @() sparse(AA_i, AA_j, ones(length(AA_i), 1), size(AA, 1), size(AA, 2));
                option.cl=[0;0;-inf(n,1)];
                option.cu=[0;0;zeros(n,1)];
            case 'DCP3'
                x0ipopt = [x0.x;x0.y;x0.w];
                funchandle.objective = @(x) obj_nlp3(x, n); % output first term only
                funchandle.gradient = @(x) con_nlp3(x, n); % output the second term only
                % Set variable bounds
                option.lb = zeros(3 * n, 1);
                option.ub = inf(3 * n, 1);

                % Set linear constraints
                Aeq = [B, -A, -eye(n);
                    ones(1, n), zeros(1, n), zeros(1, n)];
                beq = [zeros(n, 1); 1];
                funchandle.constraints=@(x) Aeq*x - beq;
                funchandle.jacobian=@(x) sparse(Aeq);
                [Aeq_i, Aeq_j] = find(Aeq);
                funchandle.jacobianstructure = @() sparse(Aeq_i, Aeq_j, ones(length(Aeq_i), 1), size(Aeq, 1), size(Aeq, 2));
                option.cl=zeros(n+1,1);
                option.cu=zeros(n+1,1);
        end
        option.ipopt.print_level           = 0;
        option.ipopt.hessian_approximation = 'limited-memory';

        tic
        [xopt,info] = ipopt(x0ipopt,funchandle,option);

        % get results
        cputime=toc;
        fopt=funchandle.objective(xopt);
        iters = info.iter;
        fprintf('Solution of IPOPT for model %s: time %.3f sec, obj %.5e iters %d\n',modelname,cputime,fopt,iters);
        c=computerr_disp(xopt,n,A,B);
end