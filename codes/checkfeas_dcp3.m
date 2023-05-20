function yesno = checkfeas_dcp3(obj,n,zk,tol_infeas)
xx=zk(1:n);
yy=zk(n+1:2*n);
ww=zk(2*n+1:3*n);
if all(xx>=-tol_infeas) && all(yy>=-tol_infeas) && all(ww>=-tol_infeas) ...
        && norm(obj.B*xx-obj.A*yy-ww)<=tol_infeas && abs(sum(xx)-1)<tol_infeas
    yesno = true;
else
    yesno = false;
end
end