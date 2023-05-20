function [dHx,dHy,dHw,dHz] = dH1(x,y,w,z)
xpy = x+y;
xmy = x-y;
xmwo2 = (x-w)/2;
zp1 = z+1;
zm1 = z-1;
zp1xpy = zp1^2 + xpy'*xpy;
zm1ymx = zm1^2 + xmy'*xmy;

dHx = (zp1xpy*xpy + zm1ymx*xmy)/4 + xmwo2 + 2*(x.'*x)*x;
dHy = (zp1xpy*xpy - zm1ymx*xmy)/4;
dHw = -xmwo2;
dHz = (zp1xpy*zp1 + zm1ymx*zm1)/4 + 2*z^3;
end