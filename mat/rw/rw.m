% Reifschneider-Williams rule
clear all;
damp = 1.0;

bet = 0.9925;
kap = 0.01;
phi = 10.0;
phi = 100.0;
phi = 200.0;
phi = 5.0;
% phi = 2.5;
rstar = (1/bet-1)*100;
zeta = 1.0;

pH = 0.0;
pL = 0.8;
sH = rstar;
sL = rstar-2.5;
%sL = rstar-1.9557;
% sH = 0;
% sL = -1.9557;

% initial guess for zeta = 0
% (semi-)analytical
A = [-1+(1-pH) pH -(phi-1)*(1-pH) -(phi-1)*pH;
    kap 0 -1+bet*(1-pH) bet*pH;
    (1-pL) -1+pL (1-pL) pL;
    0 kap bet*(1-pL) -1+bet*pL];
b = [rstar-sH;0;-sL;0];
x = A\b;

yH  = x(1);
yL  = x(2);
piH = x(3);
piL = x(4);
% check>0
iH = rstar + phi*((1-pH)*piH + pH*piL);

% policy function iteration
N = 101;
Ns = 2;
ymat0 = ones(N,1)*[yH yL];
pmat0 = ones(N,1)*[piH piL];
imat0 = ones(N,1)*[iH 0];
% ymat0 = zeros(N,2);
% pmat0 = zeros(N,2);
% imat0 = zeros(N,2);
mmat0 = zeros(N,2);
nmat0 = imat0; % nortional rate
ymat1 = ymat0;
pmat1 = pmat0;
imat1 = imat0;
mmat1 = mmat0;
nmat1 = nmat0;

mmax = 0.0;
% mmin = -2.0;
% mmin = -4.0;
mmin = -6.0;
mmin = -40.0;
mmin = -60.0;
mmin = -4.0;
% mmin = -2.0;
mgrid = linspace(mmin,mmax,N)';

Gs = [sH; sL];
Ps = [1-pH pH; 1-pL pL]; 

invlam = 1/0.003;
sig = 1.0;
[xvec0 pvec0 rvec0 vvec0] = calcmpe(Gs,Ps,bet,kap,invlam,sig,rstar);
disp([xvec0(2) pvec0(2)]);
pause;

crit = 1e-4;
diff = 1e+4;
iter = 0;
options = optimset('TolFun',1e-10);

while(diff>crit)

    for is = 1:Ns

        s0 = Gs(is);
        
        eyvec = zeros(N,1);
        epvec = zeros(N,1);
        
        for js = 1:Ns

            eyvec = eyvec + Ps(is,js)*ymat0(:,js);
            epvec = epvec + Ps(is,js)*pmat0(:,js);
            
        end
        
        for im = 1:N            

            mpast = mgrid(im);
            
            % solve for m0
            x = fzero(@rw_eqm,mpast,options,mpast,mgrid,epvec,rstar,phi,zeta);
%             x = golden('rw_eqm',mmin,mmax,1e-10,mpast,mgrid,epvec,rstar,phi,zeta);
            m0 = x;
            
            ey0 = intf1(mgrid,eyvec,m0);
            ep0 = intf1(mgrid,epvec,m0);
            is0 = rstar + phi*ep0;
            i0 = max(0, is0+mpast);
            y0 = ey0 - (i0 - ep0 - s0);
            p0 = kap*y0 + bet*ep0;
            if (is0>0.0)
                m0 = mpast - (i0-is0);
            else
                m0 = mpast - zeta*(0-is0); 
            end

            ymat1(im,is) = y0;
            pmat1(im,is) = p0;
            imat1(im,is) = i0; %i0;
            mmat1(im,is) = m0;
            nmat1(im,is) = is0+mpast; %i0;
            
        end

    end
    
    ydiff = max(max(abs(ymat1-ymat0)));
    pdiff = max(max(abs(pmat1-pmat0)));
%    idiff = max(max(abs(imat1-imat0)));
    diff = max([ydiff pdiff]);
    iter = iter + 1;
    disp([iter diff]);

    ymat0 = damp*ymat1 + (1-damp)*ymat0;
    pmat0 = damp*pmat1 + (1-damp)*pmat0;
    imat0 = damp*imat1 + (1-damp)*imat0;
    mmat0 = damp*mmat1 + (1-damp)*mmat0;
    nmat0 = damp*nmat1 + (1-damp)*nmat0;
    
end

save pf.mat Gs Ps mgrid bet kap phi rstar zeta ymat0 pmat0 imat0 mmat0 nmat0;