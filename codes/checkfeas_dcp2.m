function yesno = checkfeas_dcp2(obj,n,zk,tol_infeas)
xx=zk(1:n);
yy=zk(n+1:2*n);
zz=zk(2*n+1);
if all(xx>=-tol_infeas) && all(yy>=-tol_infeas) && all(zz>=-tol_infeas) ...
        && all(obj.B*xx-obj.A*yy>=-tol_infeas) && abs(sum(xx)-1)<tol_infeas && ...
        abs(sum(yy)-zz)<tol_infeas
    yesno = true;
else
    yesno = false;
end
end