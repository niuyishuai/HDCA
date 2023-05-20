function bestz = findbestsz1(coef,ubz)
% coef = [a1,a2,a3,a4]; ubz>0.

rr=roots(coef(1:4));
rr = real(rr(abs(imag(rr))<1e-8));
r1=rr(rr>=0);
r2=r1(r1<=ubz);
L=[0;ubz;r2(:)];
vv=polyval(coef./[4,3,2,1,1],L);
[~,idx]=min(vv);
bestz=min(L(idx));
%bestz=max(L(idx)); % for largest stepsize
end