function [dHx,dHy,dHw] = dH3(x,y,w,eta)
xyxx = (x'*y)/(x'*x);

dHx = eta*x + 2*xyxx*y - 2*xyxx^2*x + (x-w)/2;
dHy = eta*y + 2*xyxx*x;
dHw = (w-x)/2;
end