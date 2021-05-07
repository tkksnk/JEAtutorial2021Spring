function [xvec0 pvec0 rvec0 vvec0] = calcmpe(Gs,Ps,bet,kap,invlam,sig,rstar)

sh = Gs(1);
sl = Gs(2);
ph = Ps(1,2);
pl = Ps(2,2);

Ax = zeros(6,6);
Ax(1,:) = [sig*ph -(1-ph) 1 -sig*ph -ph 0];
Ax(2,:) = [-kap 1-bet*(1-ph) 0 0 -bet*ph 0];
Ax(3,:) = [-1 -kap*invlam 0 0 0 0];
Ax(4,:) = [-sig*(1-pl) -(1-pl) 0 sig*(1-pl) -pl 0];
Ax(5,:) = [0 -bet*(1-pl) 0 -kap 1-bet*pl 0];
Ax(6,:) = [0 0 0 -1 -kap*invlam -1];
bx = zeros(6,1);
bx(1) = sh-rstar;
bx(4) = sl;
xx = inv(Ax)*bx;
yh = xx(1);
yl = xx(4);
pih = xx(2);
pil = xx(5);
rh = xx(3);
rl = -rstar;

Av = zeros(2,2);
Av(1,:) = [1-bet*(1-ph) -bet*ph];
Av(2,:) = [-bet*(1-pl) 1-bet*pl];
bv = zeros(2,1);
bv(1) = -yh^2 - invlam*pih^2;
bv(2) = -yl^2 - invlam*pil^2;
vv = inv(Av)*bv;
vh = vv(1);
vl = vv(2);
xvec0 = [yh yl]';
pvec0 = [pih pil]';
rvec0 = [rh rl]';
vvec0 = [vh vl]';