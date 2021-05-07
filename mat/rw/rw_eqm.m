function f = rw_eqm(x,mpast,mgrid,epvec,rstar,phi,zeta)

m0 = x;
epi = intf1(mgrid,epvec,m0);
istar = rstar + phi*epi;
i0 = max(istar+mpast,0);

if (istar<0.0)
    f = -m0 + mpast - zeta*(0-istar); 
else
    f = -m0 + mpast - (i0-istar);
end

% f = abs(f);
