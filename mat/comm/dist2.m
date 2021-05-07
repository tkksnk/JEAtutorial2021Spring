function f = dist2(x0,sh,Ps,bet,invlam,sig,rstar,xdtar,pdtar)

sl = x0(1);
kap = x0(2);

[xvec0 pvec0] = calcmpe([sh sl]',Ps,bet,kap,invlam,sig,rstar);

f = (xvec0(2)-xdtar)^2 + (pvec0(2)-pdtar)^2;