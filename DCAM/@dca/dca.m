classdef dca < matlab.mixin.Copyable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DC programming solver
    %
    % Author: Yi-Shuai NIU
    % History:
    % - code initialized 2019-4 for general DCA
    % - modified in 2020 for BDCA
    % - modified in 2021 for InDCA
    % - modified in 2023 for ADCA, HDCA-LI and HDCA-NI for solving AEiCP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        dcp % dc program
        x0 % initial point
    end
    properties(GetAccess = public,SetAccess = private) % read only
        fopt = inf % objective value
        xopt = [] % optimal solution
        iter = 0 % iterations of dca
    end
    properties
        tolf = 1e-6 % tolerence for objective function
        tolx = 1e-6 % tolerence for iterative point
        maxiter = 1e+4 % max iterations for dca
        verbose = true  %true: display iterations of dca, false: otherwise
        convexsolver = '' % convex subproblem solver, default quadprog
        %convexsolver_option = optimset('Display','off') %
        % problem data
        A;
        B;
        model; % DCP1|DCP2|DCP3
        % parameters for plotting
        plot = true % true: draw iterations of dca, false: otherwise
        plotlinetype='b-s' % this option set the plot line type
        plotsemilogy=false % if true plot semilogy
        plotinnewfig=true % if true, plot in a new fig, otherwise, plot in fig 1
        % strong convexification parameter
        strongcvx=0.1; % set 0 to ignore this term 
        % parameters for line search
        linesearch = true % use line search for acceleration
        linesearch_type = 'exact' %'armijo|exact'
        alphabar=1; % upper bound of line search
        % parameters for armijo line search
        armijo_beta=0.1;
        armijo_sigma=1e-3;
        armijo_tol=1e-8;
        % parameters for local solver (MOSEK)
        localsol_tol=1e-8;
        % parameters for nesterov acceleration
        nesterov=false; % use nesterov acceleration
        restartperiod=inf; % period for restarting nesterov
        adca_q=10; % parameter q for adca
        % parameters for inertial acceleration
        inertial=false; % use inertial acceleration
        inertial_conserv=false; % convervative option for inertial force
        % other options
        savepoints=false; % if true then we will save objective values of all iterations
    end

    methods
        function obj = dca(dcp,x0)
            % dca constructor
            % obj = dca(dcp,x0)
            % where dcp is a dc program object
            % x0 is an initial point. It will be a random point if x0 is not given.
            if nargin==1 % if x0 is not given
                obj.dcp = dcp;
                x0 = rand(size(dcp.X));
                obj.x0 = x0(:);
            elseif nargin==2 % if x0 is given
                obj.dcp = dcp;
                obj.x0 = x0(:);
            else
                error('wrong input arguments.');
            end
        end
        function xopt = get.xopt(obj)
            % get optimal solution
            xopt = obj.xopt;
        end
        function fopt = get.fopt(obj)
            % get objective value
            fopt = obj.fopt;
        end
        function set.tolf(obj,val)
            % set tolerence of objective function
            obj.tolf = val;
        end
        function set.tolx(obj,val)
            % set tolerence of iterative point
            obj.tolx = val;
        end
        function set.maxiter(obj,val)
            % set max iterations of dca
            obj.maxiter = val;
        end
        function set.plot(obj,val)
            % set plotting option of dca
            obj.plot = val;
        end
        function set.verbose(obj,val)
            % set verbose option of dca
            obj.verbose = val;
        end
        function set.convexsolver(obj,val)
            % set convex subproblem solver (used for yalmip)
            obj.convexsolver = val;
        end
        function set.linesearch(obj,yn)
            % set line search
            obj.linesearch = yn;
        end
        % dca algorithm for solving dcp with starting point x0
        function status = optimize(obj)
            % dca optimizer
            % status.flag : 0 dca converges with tolx or tolf
            %               1 maxiter exceed
            %               2 problem infeasible or unbounded
            % status.info : solution informations.
            % status.iter : number of iterations of dca.
            % status.time : cpu time (sec.)
            % status.avgt : average time for each iteration (sec.)

            mosekcmd='minimize echo(0)';
            if obj.savepoints
                fobjlst = zeros(obj.maxiter,1);
                clst = zeros(obj.maxiter,1);
                timelst = fobjlst;
            end
            % plotting if actived
            if (obj.plot==1)
                if (obj.plotinnewfig)
                    figure
                else
                    figure(1);
                end
            end
            % display banner message
            if (obj.verbose == 1)
                fprintf('------------------------------------------------------------\n');
                fprintf('DCA-type algorithms (version 2.0) \nCopyright (c) Yi-Shuai Niu 2023\nPlateform: %s\n',computer);
                if obj.linesearch
                    switch obj.linesearch_type
                        case 'armijo'
                            fprintf('* activate Armijo linesearch acceleration\n');
                        case 'exact'
                            fprintf('* activate exact linesearch acceleration\n');
                    end
                end
                if obj.nesterov
                    fprintf('* activate Nesterov acceleration\n');
                end
                if obj.inertial
                    fprintf('* activate Inertial acceleration\n');
                end
                fprintf('------------------------------------------------------------\n');
                tabw = [5,12,12,12];
                fprintf('%-*s | %-*s | %-*s | %-*s\n',tabw(1),'iter',tabw(2),'fobj',tabw(3),'dx',tabw(4),'df');
            end
            cputime=tic;
            % different models
            switch obj.model
                case 'DCP1'
                    obj.iter = 0;
                    n = numel(obj.x0.x);
                    xk = [obj.x0.x;obj.x0.y;obj.x0.w;obj.x0.z];
                    fk=obj.dcp.F.f(xk,n,obj.A,obj.B,1); % 0 for computing both f(xk) and df(xk), 1 for computing f(xk) only

                    [subprob,param_mosek] = initsubprob1(n,obj.A,obj.B,obj.localsol_tol);
                    subprob=updatesubprob_strongconvex(subprob,obj.strongcvx*speye(3*n+1));
                    xk = [xk;...
                        (obj.x0.z+1)^2;...
                        norm(obj.x0.y-obj.x0.x)^2;...
                        (obj.x0.z-1)^2;...
                        norm(obj.x0.y+obj.x0.x)^2;...
                        obj.x0.z^2;...
                        obj.x0.x'*obj.x0.x
                        ];
                    dk=0;
                    if obj.inertial && ~obj.nesterov
                        % set stepsize for inertial force. Different setting will be used to nesterov
                        gamma = obj.strongcvx; 
                    else
                        gamma = 0;
                    end

                    % Main loop
                    while obj.iter < obj.maxiter
                        xkorg=xk; % avoid modifications by Nesterov
                        if obj.nesterov && obj.iter<1
                            % initialize nesterov parameters
                            tk=1;
                            fqlst=-inf(obj.adca_q,1);
                            q=1;
                        elseif obj.nesterov 
                            if mod(obj.iter,obj.restartperiod)==0
                                % restart nesterov parameter
                                tk=1;
                            end
                            % Nesterov's acceleration
                            tk1 = (1+sqrt(1+4*tk^2))/2;
                            betak = (tk-1)/tk1;
                            gamma = 2*obj.strongcvx*(1-betak^2)/(3-betak^2);
                            zk = xk + ((tk-1)/(tk1))*dk;
                            % ADCA only
                            if ~obj.inertial
                                if obj.dcp.F.f(zk,n,obj.A,obj.B,0) <= max(fqlst)
                                    xk = zk;
                                end
                            else % HDCA-NI only
                                if obj.dcp.F.f(zk,n,obj.A,obj.B,0) + (2*obj.strongcvx-gamma)*normdk^2/4 <= max(fqlst) %&& checkfeas_dcp1(obj,n,zk,1e-8)
                                    xk = zk;
                                end
                            end
                            tk = tk1;
                        end

                        % update subgradient
                        [dHx,dHy,dHw,dHz] = dH1(xk(1:n),xk(n+1:2*n),xk(2*n+1:3*n),xk(3*n+1)); % compute dH1(xk)
                        % initialize mosek problem
                        if (obj.inertial && obj.inertial_conserv && obj.iter<2) || (~obj.inertial)
                            dd=obj.strongcvx*xk;
                        else
                            dd=obj.strongcvx*xk+gamma*dk;
                        end
                        % the linear part of QP formulation without u.
                        subprob.c=-[[dHx;dHy;dHw;dHz]+dd(1:3*n+1);zeros(6,1)];
                        subprob.barx0 = xk;

                        % solve convex subproblem via mosek
                        [~,res] = mosekopt(mosekcmd,subprob,param_mosek);
                        xk1=res.sol.itr.xx;
                        obj.iter = obj.iter+1;
                        %res.sol.itr.pobjval
                        %[xk1,fkval1] = solvecq1_cvx(dHx+dd(1:n),dHy+dd(n+1:2*n),dHw+dd(2*n+1:3*n),dHz+dd(3*n+1),obj.A,obj.B,obj.strongcvx);
                        %fkval1
                        %res_msk=solvecq1_msk(dHx+dd(1:n),dHy+dd(n+1:2*n),dHw+dd(2*n+1:3*n),dHz+dd(3*n+1),obj.A,obj.B,obj.strongcvx);
                        
                        [fk1,dfk1] = obj.dcp.F.f(xk1,n,obj.A,obj.B,0);
                        dk = xk1 - xkorg;
                        normdk = norm(dk);

                        if obj.nesterov
                            if mod(q,obj.adca_q+1)==0
                                q=1;
                            end
                            % ADCA only
                            if ~obj.inertial
                                fqlst(q)=fk1;
                            else % HDCA-NI
                                fqlst(q)=fk1 + (2*obj.strongcvx-gamma)*normdk^2/4;
                            end
                            q=q+1;
                        end

                        % accelerate with line search
                        if obj.linesearch == true
                            dk3n = dk(1:3*n+1);
                            if (dk3n'*dfk1 < 0 )

                                %activesetxk = checkactiveset(xk);
                                %activesetxk1 = checkactiveset(xk1);
                                %if (min(activesetxk - activesetxk1)>=0 && dk3n'*dfk1 < 0 )
                                switch obj.linesearch_type
                                    case 'armijo'
                                        [xacc,facc] = armijo(obj,n,xk1,fk1,dk,obj.armijo_beta,obj.armijo_sigma,obj.alphabar,obj.armijo_tol);
                                        %[xacc,facc,obj.stepsize]=armijo_adaptive(obj,fk1,xk1,dk,obj.stepsize); % increase many times NEED TO CHECK
                                        if facc < fk1
                                            if obj.verbose == 1
                                                fprintf('accelerated: reduced %17.3e  moved %17.3e \n',facc-fk1,norm(xacc-xk1));
                                            end
                                            xk1 = xacc;
                                            fk1 = facc;
                                        end
                                    case 'exact'
                                        Ik = dk3n<0;
                                        dkx = dk(1:n);
                                        dky = dk(n+1:2*n);
                                        dkw = dk(2*n+1:3*n);
                                        dkz = dk(3*n+1);
                                        vkx = xk1(1:n);
                                        vky = xk1(n+1:2*n);
                                        vkw = xk1(2*n+1:3*n);
                                        vkz = xk1(3*n+1);
                                        sig0=dkz*dkx;
                                        sig1=dky-dkz*vkx-vkz*dkx;
                                        sig2=vky-vkz*vkx;
                                        ploycoefs = zeros(1,5);
                                        ploycoefs(1)=4*dkz^2*norm(dkx)^2; %a1
                                        ploycoefs(2)=-6*sig0'*sig1; %a2
                                        ploycoefs(3)=2*(norm(sig1)^2 + dkw'*dkx - 2*sig0'*sig2); %a3
                                        ploycoefs(4)=2*sig2'*sig1 + vkw'*dkx + vkx'*dkw; %a4
                                        ploycoefs(5)=vkx'*vkw + norm(sig2)^2; %a5
                                        if sum(Ik)==0 % Ik is empty
                                            alphak = findbestsz1(ploycoefs,inf);
                                        else
                                            alphak = findbestsz1(ploycoefs,min(abs(xk1(Ik)./dk(Ik))));
                                        end
 
                                        if alphak > 0 % if alphak is not too small, update xk1
                                            xacc = xk1 + alphak*dk;
                                            facc = obj.dcp.F.f(xacc,n,obj.A,obj.B,1);
                                            if facc < fk1
                                                if obj.verbose == 1
                                                    fprintf('* acceleration: f-reduced %.3e | x-moved %.3e \n',facc-fk1,norm(xacc-xk1));
                                                end
                                                xk1 = xacc;
                                                fk1 = facc;
                                            end
                                        end

                                end

                            end
                        end

                        % compute errors
                        normx = normdk;
                        normf = abs(fk1-fk);
                        % display iterations of dca if actived
                        if (obj.verbose == 1)
                            fprintf('%-*d | %-*.5e | %-*.5e | %-*.5e\n',tabw(1),obj.iter,tabw(2),fk1,tabw(3),normx,tabw(4),normf);
                        end
                        % plotting if actived
                        if (obj.plot)
                            myplotf(fk,fk1,obj.iter,obj.plotlinetype,obj.plotsemilogy);
                        end
                        % save history
                        if obj.savepoints
                            fobjlst(obj.iter)=fk1;
                            clst(obj.iter)=computerr(xk,n,obj.A,obj.B);
                            timelst(obj.iter)=toc(cputime);
                        end
                        % check stopping
                        if (normx < obj.tolx*(1+norm(xk1)) || normf < obj.tolf*(1+abs(fk1)))
                            if (obj.verbose == 1)
                                fprintf('------------------------------------------------------------\n');
                            end
                            obj.fopt = fk1;
                            obj.xopt = xk1;
                            status = setstatus(toc(cputime),obj.iter,0,'Successfully solved.');
                            if obj.savepoints
                                status.fvallst=fobjlst(1:obj.iter);
                                status.clst=clst(1:obj.iter);                                
                                status.timelst=timelst(1:obj.iter);
                            end
                            return;
                        end

                        xk = xk1;
                        fk = fk1;
                    end % end for while
                case 'DCP2'
                    obj.iter = 0;
                    n = numel(obj.x0.x);
                    xk = [obj.x0.x;obj.x0.y;obj.x0.z];
                    fk=obj.dcp.F.f(xk,n,obj.A,obj.B,1); % 0 for computing both f(xk) and df(xk), 1 for computing f(xk) only

                    [subprob,param_mosek] = initsubprob2(n,obj.A,obj.B,obj.localsol_tol);
                    %if obj.inertial
                    subprob=updatesubprob_strongconvex(subprob,obj.strongcvx*speye(2*n+1));
                    xk = [xk;...
                        (obj.x0.z+1)^2;... % u1
                        norm(obj.x0.y-obj.x0.x)^2;... %u2
                        (obj.x0.z-1)^2;... %u3
                        norm(obj.x0.y+obj.x0.x)^2;... %u4
                        obj.x0.z^2;... %u5
                        obj.x0.x'*obj.x0.x %u6
                        ];
                    dk=0;
                    if obj.inertial && ~obj.nesterov
                        gamma = obj.strongcvx-0.001; % stepsize for inertial force
                    else
                        gamma = 0;
                    end
                    %end
                    while obj.iter < obj.maxiter
                        xkorg=xk; % avoid modifications by Nesterov
                        if obj.nesterov && obj.iter<1
                            % initialize nesterov parameters
                            tk=1;
                            fqlst=-inf(obj.adca_q,1);
                            q=1;
                        elseif obj.nesterov %&& obj.iter>=1
                            if mod(obj.iter,obj.restartperiod)==0
                                % restart nesterov parameter
                                tk=1;
                            end
                            % Nesterov's acceleration
                            tk1 = (1+sqrt(1+4*tk^2))/2;
                            betak = (tk-1)/tk1;
                            gamma = 2*obj.strongcvx*(1-betak^2)/(3-betak^2);
                            zk = xk + ((tk-1)/(tk1))*dk;
                            % ADCA only
                            if ~obj.inertial
                                if obj.dcp.F.f(zk,n,obj.A,obj.B,0) <= max(fqlst)
                                    xk = zk;
                                end
                            else % HDCA-NI only
                                if obj.dcp.F.f(zk,n,obj.A,obj.B,0) + (2*obj.strongcvx-gamma)*normdk^2/4 <= max(fqlst) %&& checkfeas_dcp1(obj,n,zk,1e-8)
                                    xk = zk;
                                end
                            end
                            tk = tk1;
                        end
                       
                        [dHx,dHy,dHz] = dH2(xk(1:n),xk(n+1:2*n),xk(2*n+1),obj.A); % compute dH2(xk)
                        % initialize mosek problem
                        if (obj.inertial && obj.inertial_conserv && obj.iter<2) || (~obj.inertial)
                            dd=obj.strongcvx*xk;
                        else
                            dd=obj.strongcvx*xk+gamma*dk;
                        end
                        % the linear part of QP formulation without u.
                        subprob.c=-[[dHx;dHy;dHz]+dd(1:2*n+1);zeros(6,1)];
                        subprob.barx0 = xk;

                        % solve convex subproblem via mosek
                        [~,res] = mosekopt(mosekcmd,subprob,param_mosek);
                        %res.sol.itr.pobjval
                        %[xk1,fkval1] = solvecq2_cvx(dHx+dd(1:n),dHy+dd(n+1:2*n),dHz+dd(2*n+1),obj.A,obj.B,obj.strongcvx);
                        %fkval1
                        if res.rcode~=0
                            fprintf('子问题没有正常求解, iter: %d, rcode: %d, %s\n',obj.iter,res.rcode,res.rcodestr);
                            %[xk1,fkval1] = solvecq1_cvx(dHx+dd(1:n),dHy+dd(n+1:2*n),dHw+dd(2*n+1:3*n),dHz+dd(3*n+1),obj.A,obj.B,obj.strongcvx);
                        end
                        xk1=res.sol.itr.xx;
                        obj.iter = obj.iter+1;

                        [fk1,dfk1] = obj.dcp.F.f(xk1,n,obj.A,obj.B,0);
                        dk = xk1 - xkorg;
                        normdk = norm(dk);

                        if obj.nesterov
                            if mod(q,obj.adca_q+1)==0
                                q=1;
                            end
                            % ADCA only
                            if ~obj.inertial
                                fqlst(q)=fk1;
                            else
                                fqlst(q)=fk1 + (2*obj.strongcvx-gamma)*normdk^2/4;
                            end
                            q=q+1;
                        end



                        % accelerate with line search
                        if obj.linesearch == true
                            dk3n = dk(1:2*n+1);
                            if (dk3n'*dfk1 < 0 )

                                %activesetxk = checkactiveset(xk);
                                %activesetxk1 = checkactiveset(xk1);
                                %if (min(activesetxk - activesetxk1)>=0 && dk3n'*dfk1 < 0 )
                                switch obj.linesearch_type
                                    case 'armijo' % TODO
                                        [xacc,facc] = armijo(obj,n,xk1,fk1,dk,obj.armijo_beta,obj.armijo_sigma,obj.alphabar,obj.armijo_tol);
                                        %[xacc,facc,obj.stepsize]=armijo_adaptive(obj,fk1,xk1,dk,obj.stepsize); % increase many times NEED TO CHECK
                                        if facc < fk1
                                            if obj.verbose == 1
                                                fprintf('accelerated: reduced %17.3e  moved %17.3e \n',facc-fk1,norm(xacc-xk1));
                                            end
                                            xk1 = xacc;
                                            fk1 = facc;
                                        end
                                    case 'exact'
                                        dkx = dk(1:n);
                                        dky = dk(n+1:2*n);
                                        dkz = dk(2*n+1);
                                        vkx = xk1(1:n);
                                        vky = xk1(n+1:2*n);
                                        vkw = obj.B*vkx - obj.A*vky;
                                        vkz = xk1(2*n+1);
                                        sig0=dkz*dkx;
                                        sig1=dky-dkz*vkx-vkz*dkx;
                                        sig2=obj.B*dkx - obj.A*dky;
                                        sig3=vky-vkz*vkx;
                                        Ik = dk3n < 0;
                                        Jk = sig2 < 0;
                                        ploycoefs = zeros(1,5);
                                        ploycoefs(1)=4*(sig0'*sig0); %a1
                                        ploycoefs(2)=-6*sig0'*sig1; %a2
                                        ploycoefs(3)=2*((sig1'*sig1) + sig2'*dkx - 2*sig0'*sig3); %a3
                                        ploycoefs(4)=2*sig3'*sig1 + vkw'*dkx + vkx'*sig2; %a4
                                        ploycoefs(5)=vkx'*vkw + norm(sig3)^2; %a5
                                        if sum(Ik)==0 % Ik is empty
                                            alphak = findbestsz1(ploycoefs,inf);
                                        else
                                            alphak = findbestsz1(ploycoefs,min([-xk1(Ik)./dk(Ik);-vkw(Jk)./sig2(Jk)]));
                                        end

                                        if alphak > 1e-8 % if alphak is not too small, update xk1
                                            xacc = xk1 + alphak*dk;
                                            facc = obj.dcp.F.f(xacc,n,obj.A,obj.B,1);
                                            if facc < fk1
                                                if obj.verbose == 1
                                                    fprintf('* acceleration: f-reduced %.3e | x-moved %.3e \n',facc-fk1,norm(xacc-xk1));
                                                end
                                                xk1 = xacc;
                                                fk1 = facc;
                                            end
                                        end

                                end

                            end
                        end

                        % compute errors
                        normx = normdk;
                        normf = abs(fk1-fk);
                        % display iterations of dca if actived
                        if (obj.verbose == 1)
                            fprintf('%-*d | %-*.5e | %-*.5e | %-*.5e\n',tabw(1),obj.iter,tabw(2),fk1,tabw(3),normx,tabw(4),normf);
                        end
                        % plotting if actived
                        if (obj.plot)
                            myplotf(fk,fk1,obj.iter,obj.plotlinetype,obj.plotsemilogy);
                        end
                        % save history
                        if obj.savepoints
                            fobjlst(obj.iter)=fk1;
                            clst(obj.iter)=computerr(xk,n,obj.A,obj.B);
                            timelst(obj.iter)=toc(cputime);
                        end
                        % check stopping
                        if (normx < obj.tolx*(1+norm(xk1)) || normf < obj.tolf*(1+abs(fk1)))
                            if (obj.verbose == 1)
                                fprintf('------------------------------------------------------------\n');
                            end
                            obj.fopt = fk1;
                            obj.xopt = xk1;
                            status = setstatus(toc(cputime),obj.iter,0,'Successfully solved.');
                            if obj.savepoints
                                status.fvallst=fobjlst(1:obj.iter);
                                status.clst=clst(1:obj.iter);                                
                                status.timelst=timelst(1:obj.iter);
                            end
                            return;
                        end

                        xk = xk1;
                        fk = fk1;
                    end % end for while
                case 'DCP3'
                    obj.iter = 0;
                    n = numel(obj.x0.x);
                    xk = [obj.x0.x;obj.x0.y;obj.x0.w];
                    fk=obj.dcp.F.f(xk,n,obj.A,obj.B,1); % 0 for computing both f(xk) and df(xk), 1 for computing f(xk) only

                    [subprob,param_mosek] = initsubprob3(n,obj.A,obj.B,obj.localsol_tol);
                    % solve linear problem via mosek to get M
                    [~,res] = mosekopt(mosekcmd,subprob,param_mosek);
                    eta = 3.2 + 20*n*res.sol.itr.pobjval^2;
                    % set the quadratic part Qo (lower triangular only) of G3
                    matz = sparse(3*n,3*n);
                    matid = speye(n);
                    Qo=matz;
                    Qo(1:n,1:n)=(eta + 1/2)*matid;
                    Qo(2*n+1:3*n,1:n)=matid/2;
                    Qo(n+1:2*n,n+1:2*n)=(2 + eta)*matid;
                    Qo(2*n+1:3*n,2*n+1:3*n)=matid/2;
                    [subprob.qosubi,subprob.qosubj,subprob.qoval] = find(tril(Qo));

                    % convert to strongly convex
                    subprob=updatesubprob_strongconvex(subprob,obj.strongcvx*speye(3*n));
                    dk=0;
                    % initialize qpsubp
                    %QP = initqpsubp3(n,eta,obj.strongcvx,obj.A,obj.B);
                    if obj.inertial && obj.linesearch % HDCA-LI
                        gamma = 2*obj.strongcvx/(1+(1+obj.alphabar)^2);
                    elseif obj.inertial && ~obj.linesearch %~obj.nesterov
                        gamma = obj.strongcvx; % stepsize for inertial force
                    else
                        gamma = 0;
                    end
                    
                    while obj.iter < obj.maxiter
                        xkorg=xk; % avoid modifications by Nesterov
                        if obj.nesterov && obj.iter<1
                            % initialize nesterov parameters
                            tk=1;
                            fqlst=-inf(obj.adca_q,1);
                            q=1;
                        elseif obj.nesterov %&& obj.iter>=1
                            if mod(obj.iter,obj.restartperiod)==0
                                % restart nesterov parameter
                                tk=1;
                            end
                            % Nesterov's acceleration
                            tk1 = (1+sqrt(1+4*tk^2))/2;
                            betak = (tk-1)/tk1;
                            gamma = 2*obj.strongcvx*(1-betak^2)/(3-betak^2);
                            zk = xk + ((tk-1)/(tk1))*dk;
                            % ADCA only
                            if ~obj.inertial
                                if obj.dcp.F.f(zk,n,obj.A,obj.B,0) <= max(fqlst)
                                    xk = zk;
                                end
                            else % HDCA-NI only
                                if obj.dcp.F.f(zk,n,obj.A,obj.B,0) + (2*obj.strongcvx-gamma)*normdk^2/4 <= max(fqlst) %&& checkfeas_dcp1(obj,n,zk,1e-8)
                                    xk = zk;
                                end
                            end
                            tk = tk1;
                        end

                        [dHx,dHy,dHw] = dH3(xk(1:n),xk(n+1:2*n),xk(2*n+1:3*n),eta); % compute dH2(xk)
                        % initialize mosek problem
                        if (obj.inertial && obj.inertial_conserv && obj.iter<2) || (~obj.inertial)
                            dd=obj.strongcvx*xk;
                        else
                            dd=obj.strongcvx*xk+gamma*dk;
                        end
                        % the linear part of QP formulation without u.
                        subprob.c=-[dHx;dHy;dHw]-dd;
                        %subprob.barx0 = xk;

                        switch obj.convexsolver
                            case 'qp'
                                % solve qp subp
                                %[xk1,fval1,niter]=qpsubp(xk,subprob.c,QP.Q,QP.Aeq,QP.beq,QP.Aineq,QP.bineq,obj.localsol_tol);
                            otherwise
                                % solve convex subproblem via mosek
                                [~,res] = mosekopt(mosekcmd,subprob,param_mosek);
                                if res.rcode~=0
                                    fprintf('子问题没有正常求解, iter: %d, rcode: %d, %s\n',obj.iter,res.rcode,res.rcodestr);
                                    %[xk1,fkval1] = solvecq1_cvx(dHx+dd(1:n),dHy+dd(n+1:2*n),dHw+dd(2*n+1:3*n),dHz+dd(3*n+1),obj.A,obj.B,obj.strongcvx);
                                end
                                xk1=res.sol.itr.xx;
                        end
                        obj.iter = obj.iter+1;
                        
                        [fk1,dfk1] = obj.dcp.F.f(xk1,n,obj.A,obj.B,0);
                        dk = xk1 - xkorg;
                        normdk = norm(dk);

                        if obj.nesterov
                            if mod(q,obj.adca_q+1)==0
                                q=1;
                            end
                            % ADCA only
                            if ~obj.inertial
                                fqlst(q)=fk1;
                            else
                                fqlst(q)=fk1 + (2*obj.strongcvx-gamma)*normdk^2/4;
                            end
                            q=q+1;
                        end

                        % accelerate with line search
                        if obj.linesearch == true
                            if (dk'*dfk1 < 0 )

                                %activesetxk = checkactiveset(xk);
                                %activesetxk1 = checkactiveset(xk1);
                                %if (min(activesetxk - activesetxk1)>=0 && dk3n'*dfk1 < 0 )
                                [xacc,facc] = armijo(obj,n,xk1,fk1,dk,obj.armijo_beta,obj.armijo_sigma,obj.alphabar,obj.armijo_tol);
                                %[xacc,facc,obj.stepsize]=armijo_adaptive(obj,fk1,xk1,dk,obj.stepsize); % increase many times NEED TO CHECK
                                if facc < fk1
                                    if obj.verbose == 1
                                        fprintf('accelerated: f-reduced %.3e  | x-moved %.3e \n',facc-fk1,norm(xacc-xk1));
                                    end
                                    xk1 = xacc;
                                    fk1 = facc;
                                end
                            end
                        end

                        % compute errors
                        normx = normdk;
                        normf = abs(fk1-fk);
                        % display iterations of dca if actived
                        if (obj.verbose == 1)
                            fprintf('%-*d | %-*.5e | %-*.5e | %-*.5e\n',tabw(1),obj.iter,tabw(2),fk1,tabw(3),normx,tabw(4),normf);
                        end
                        % plotting if actived
                        if (obj.plot)
                            myplotf(fk,fk1,obj.iter,obj.plotlinetype,obj.plotsemilogy);
                        end
                        % save history
                        if obj.savepoints
                            fobjlst(obj.iter)=fk1;
                            clst(obj.iter)=computerr(xk,n,obj.A,obj.B);
                            timelst(obj.iter)=toc(cputime);
                        end
                        % check stopping
                        if (normx < obj.tolx*(1+norm(xk1)) || normf < obj.tolf*(1+abs(fk1)))
                            if (obj.verbose == 1)
                                fprintf('------------------------------------------------------------\n');
                            end
                            obj.fopt = fk1;
                            obj.xopt = xk1;
                            status = setstatus(toc(cputime),obj.iter,0,'Successfully solved.');
                            if obj.savepoints
                                status.fvallst=fobjlst(1:obj.iter);
                                status.clst=clst(1:obj.iter);                                
                                status.timelst=timelst(1:obj.iter);
                            end
                            return;
                        end

                        xk = xk1;
                        fk = fk1;
                    end % end for while

            end % end for switch

            % maxiter exceed
            if (obj.verbose == 1)
                fprintf('------------------------------------------------------------\n');
            end
            obj.fopt = fk1;
            obj.xopt = xk1;
            status = setstatus(toc(cputime),obj.iter,1,'Max interation exceed.');
            if obj.savepoints
                status.fvallst=fobjlst(1:obj.iter);
                status.clst=clst(1:obj.iter);    
                status.timelst=timelst(1:obj.iter);
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myplotf(fk,fk1,iter,plotline,plotsemilogy)
if iter == 1
    if plotsemilogy
        semilogy(nan,nan);
    end
    hold on;
    title('Iterations');
    xlabel('iter');
    ylabel('fval');
end
if iter>1
    %h = animatedline(gca,'Marker','s','LineStyle','-','Color','b');
    %addpoints(h,iter,fk1);
    plot([iter-1,iter], [fk,fk1],plotline);
    drawnow limitrate
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting dca solution status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = setstatus(timer,iter,flag,info)
status.time = timer;
status.iter = iter;
status.avgt = timer/iter;
status.flag = flag;
status.info = info;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% armijo
% classical armijo rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z,fz] = armijo(obj,n,v,fv,d,beta,sigma,alphabar,eps)
nd = norm(d);
idx= d<0;
alphabark=min(-v(idx)./d(idx));
alpha = min([alphabar,alphabark]);
while (alpha*nd>eps)
    z=v+alpha*d;
    if z>=-1e-8 %checkfeas(z)
        fz=obj.dcp.F.f(z,n,obj.A,obj.B,1);
        delta=fv - fz - sigma*alpha^2*nd^2;
        if delta >0
            %fprintf('alpha:%.6f\n',alpha);
            return;
        end
    end
    alpha=beta*alpha;
end
z=v;
fz=fv;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update subproblem structure by introducing strongly convex regularizer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subprob=updatesubprob_strongconvex(subprob,addQo)
[subi,subj,val]=find(tril(addQo));
subprob.qosubi=[subprob.qosubi;subi];
subprob.qosubj=[subprob.qosubj;subj];
subprob.qoval=[subprob.qoval;val];
end
