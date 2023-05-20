function [dHx,dHy,dHz] = dH2(x,y,z,A)
xpy = x+y;
xmy = x-y;
zp1 = z+1;
zm1 = z-1;
zp1xpy = zp1^2 + xpy'*xpy;
zm1ymx = zm1^2 + xmy'*xmy;
xpAy = x + A*y;

dHx = (zp1xpy*xpy + zm1ymx*xmy)/4 + 2*(x.'*x)*x + xpAy/2;
dHy = (zp1xpy*xpy - zm1ymx*xmy)/4 + (A'*xpAy)/2;
dHz = (zp1xpy*zp1 + zm1ymx*zm1)/4 + 2*z^3;
end