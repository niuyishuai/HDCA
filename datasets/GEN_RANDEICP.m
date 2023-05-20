%% generate RANDEICP

randbnd=[-1,1];
for n=[10,100,500]
    for nprobs=1:10
        T=rand(n)*(randbnd(2)-randbnd(1))+randbnd(1);
        mu = -min([0,eigs(T+T',1,'smallestreal')]);
        A=T+mu*eye(n);
        B=eye(n);

        x0.x = rand(n,1);
        x0.x = x0.x/sum(x0.x);
        x0.y = rand(n,1);
        x0.w = B*x0.x - A*x0.y;
        x0.z = sum(x0.y);

        filename = sprintf("RAND1\\RAND(1,%d)_%d.mat",n,nprobs);
        save(filename,'n','A','B','x0');

        B = 10*B;
        for i=1:n
            for j=i+1:min([i+4,n])
                B(i,j)=-1;
            end
            for j=max([1,i-4]):i-1
                B(i,j)=-1;
            end
        end

        x0.x = rand(n,1);
        x0.x = x0.x/sum(x0.x);
        x0.y = rand(n,1);
        x0.w = B*x0.x - A*x0.y;
        x0.z = sum(x0.y);

        filename = sprintf("RAND2\\RAND(2,%d)_%d.mat",n,nprobs);
        save(filename,'n','A','B','x0');
    end
end