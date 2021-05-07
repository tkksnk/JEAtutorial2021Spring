% function [mmax mmin nmax nmin sl kap] = mainsp_LOOP(xdtar,pdtar,ph,pl,invlam,mmin,mmax,nmin,nmax,pcf,N)

xdtar = -7.0;
pdtar = -1.0/4;
ph = 0.0;
pl = 0.8;
invlam = 1/0.003;

% checkMult = 1;
mmax = 0.0;
mmin = -200.0;
nmax = 50.0;
nmin = 0.0;

damp = 1.0;
% plotgraph = 1;
% prevload = 0;
diagnum = 100;
% initMPE = 0;
critc = 1e-5;
% crits = 1e-4; % >= critc

bet = 0.9925;
sig = 1.0;
%invlam = 200; %100; %6.25;  % it is necessary for invlam to be low so that \phi_1 takes neg. values.
rstar = (1/bet-1)*100;
sh = 0.0;

sl = -2.5;
kap = 0.01; %0.25/7;
% calc sl, given the target of xd
Ps = [1-ph ph; 1-pl pl];
% [x fval flag] = fminsearch(@dist2,[sl kap]',[],sh,Ps,bet,invlam,sig,rstar,xdtar,pdtar);
% sl = x(1)
% kap = x(2)
% flag
% %sl = calcsl(xdtar,sh,Ps,bet,kap,invlam,sig,rstar);
% disp(' ');
% disp(sprintf(' (invlam, ph, pl)=(%8.5f, %8.5f, %8.5f)',invlam,ph,pl));

%% set up grid points
ns = 2;
nm = 21;
nn = 21;

Gs = [sh sl]';
%Ps = [1-ph ph; 1-pl pl];

knotsm = linspace(mmin,mmax,nm)';
knotsn = logspace(log(1.0)/log(10.0), log(nmax-nmin+1.0)/log(10.0), nn)';
knotsn = knotsn + (nmin - 1.0);
% knotsn = linspace(0,nmax,nn)';

disp(' ');
disp(' computing discretionary policy...');
disp(' ');
[xvec0 pvec0 rvec0 vvec0] = calcmpe(Gs,Ps,bet,kap,invlam,sig,rstar);
eval(sprintf('save wdmat_invlam%6.5fph%1.5fpl%1.2fkap%1.5f.mat vvec0 xvec0 pvec0 rvec0 Gs Ps ns rstar',invlam,ph,pl,kap));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% commitment policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial values
nmat0 = zeros(ns,nm,nn);
mmat0 = zeros(ns,nm,nn);
xmat0 = zeros(ns,nm,nn);
pmat0 = zeros(ns,nm,nn);
rmat0 = zeros(ns,nm,nn);
vmat0 = zeros(ns,nm,nn);

crit = 1e-5;
diffin = 1e+4;
iterin = 0;

disp(' ');
disp(' computing commitment policy...');
disp(' ');
disp('          norm of  phi1     phi2     x        pai      r        v');

tic;

while (diffin>crit)

    [mmat1 nmat1 xmat1 pmat1 rmat1 vmat1] = comm_iter(Gs,Ps,knotsm,knotsn,mmat0,nmat0,xmat0,pmat0,vmat0,bet,kap,invlam,sig,rstar);

    % update
    diffm = max(max(max(abs(mmat1-mmat0))));
    diffn = max(max(max(abs(nmat1-nmat0))));
    diffx = max(max(max(abs(xmat1-xmat0))));
    diffp = max(max(max(abs(pmat1-pmat0))));
    diffr = max(max(max(abs(rmat1-rmat0))));
    diffv = max(max(max(abs(vmat1-vmat0))));
    diffin = max([diffm diffn diffx diffp diffr diffv]);
    iterin = iterin+1;
    mmat0 = mmat1;
    nmat0 = nmat1;
    xmat0 = xmat1;
    pmat0 = pmat1;
    rmat0 = rmat1;
    vmat0 = vmat1;
    
    % diagnosis
    if (mod(iterin,1)==0)
        s = sprintf( ' iteration %4d:  (%6.5f, %6.5f, %6.5f, %6.5f, %6.5f, %6.5f)', ...
            iterin, diffm, diffn, diffx, diffp, diffr, diffv);    
        disp(s);
    end
    
end

toc;
disp(' ');

vmatc = vmat0;
mmatc = mmat0;
nmatc = nmat0;
xmatc = xmat0;
pmatc = pmat0;
rmatc = rmat0;

eval(sprintf('save wcmat_invlam%6.5fph%1.5fpl%1.2fkap%1.5f.mat vmat0 mmat0 nmat0 xmat0 pmat0 rmat0 Gs Ps knotsm knotsn ns nm nn rstar',invlam,ph,pl,kap));
% % pause;
% 
% %% calculate N-period punishment values
% % NOTE: When N=1, should it be the value under MPE?
% indm = find(abs(knotsm)==min(abs(knotsm)));
% indn = find(abs(knotsn)==min(abs(knotsn)));
% vvec0 = vmat0(:,indm,indn);
% pvec0 = pmat0(:,indm,indn);
% xvec0 = xmat0(:,indm,indn);
% wvec0 = calcwdn(vvec0,pvec0,xvec0,Gs,Ps,N,bet,kap,invlam,sig,rstar);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% sustainable policy
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% initial values
% % initial guesses with MPE
% if (initMPE==1)
%     vmat0(1,:,:) = vvec0(1)*ones(1,nm,nn);
%     vmat0(2,:,:) = vvec0(2)*ones(1,nm,nn);
%     xmat0(1,:,:) = xvec0(1)*ones(1,nm,nn);
%     xmat0(2,:,:) = xvec0(2)*ones(1,nm,nn);
%     pmat0(1,:,:) = pvec0(1)*ones(1,nm,nn);
%     pmat0(2,:,:) = pvec0(2)*ones(1,nm,nn);
% end
% 
% % if (prevload==1)
% %     eval(sprintf('load wsmat_invlam%6.5fph%1.5fpl%1.2f.mat vmat0 mmat0 nmat0 xmat0 pmat0 rmat0 zmat0 idmat Gs Ps knotsm knotsn ns nm nn rstar',invlam,ph1,pl));
% % end
% 
% crit = 1e-5;
% diffout = 1e+4;
% iterout = 0;
% 
% h = figure;
% disp(' ');
% disp(' computing sustainable policy...');
% 
% tic;
% 
% while (diffout>crit)
% 
%     disp(sprintf(' punishment values: (%10.5f, %10.5f)',wvec0(1),wvec0(2)));
%     disp(' ');
%     disp('          norm of  phi1     phi2     x        pai      r        v');
% 
%     % the sustainability constraint never binds at the low state...
%     %mmat0(1,:,:) = min(0,mmat0(1,:,:));
%     vmat0 = vmatc;
%     mmat0 = mmatc;
%     nmat0 = nmatc;
%     xmat0 = xmatc;
%     pmat0 = pmatc;
%     vmat0(2,:,:) = zeros(nm,nn);
% %    vmat0(2,:,:) = zeros(nm,nn);
%     zmat0 = ones(ns,nm,nn);
%     
% %    if (iterout>0); diagnum = 1; end;
% 
%     diffin = 1e+4;
%     iterin = 0;
% 
%     while (diffin>crit)
%     
%         [mmat1 nmat1 xmat1 pmat1 rmat1 vmat1 zmat1 idmat] = ...
%             sus_iter(Gs,Ps,knotsm,knotsn,mmat0,nmat0,xmat0,pmat0,vmat0,wvec0,bet,kap,invlam,sig,rstar,checkMult);    
%         diffv = max(max(max(abs(vmat1-vmat0))));
%         diffm = max(max(max(abs(mmat1-mmat0))));
%         diffn = max(max(max(abs(nmat1-nmat0))));
%         diffx = max(max(max(abs(xmat1-xmat0))));
%         diffp = max(max(max(abs(pmat1-pmat0))));
%         diffr = max(max(max(abs(rmat1-rmat0))));
%         diffin = max([diffm diffn diffx diffp diffr diffv]);
%         iterin = iterin+1;
%         mmat0 = damp*mmat1 + (1-damp)*mmat0;
%         nmat0 = damp*nmat1 + (1-damp)*nmat0;
%         xmat0 = damp*xmat1 + (1-damp)*xmat0;
%         pmat0 = damp*pmat1 + (1-damp)*pmat0;
%         rmat0 = damp*rmat1 + (1-damp)*rmat0;
%         vmat0 = damp*vmat1 + (1-damp)*vmat0;
%         zmat0 = damp*zmat1 + (1-damp)*zmat0;
% 
%         % diagnosis
%         if (mod(iterin,diagnum)==0)
%             s = sprintf( ' iteration %4d:  (%6.5f, %6.5f, %6.5f, %6.5f, %6.5f, %6.5f)', ...
%                 iterin, diffm, diffn, diffx, diffp, diffr, diffv);    
%             disp(s);
% 
%             if (plotgraph==1)
% 
%                 figure(h);
%                 subplot(221);
%                 mesh(knotsm,knotsn,reshape(rmat0(1,:,:),[nm nn])');
%                 xlim([knotsm(1) knotsm(end)]); ylim([knotsn(1) knotsn(end)]);
%                 subplot(222);
%                 mesh(knotsm,knotsn,reshape(rmat0(ns,:,:),[nm nn])');
%                 xlim([knotsm(1) knotsm(end)]); ylim([knotsn(1) knotsn(end)]);
%                 subplot(223);
%                 mesh(knotsm,knotsn,reshape(zmat0(1,:,:),[nm nn])');
%                 xlim([knotsm(1) knotsm(end)]); ylim([knotsn(1) knotsn(end)]);
%                 subplot(224);
%                 mesh(knotsm,knotsn,reshape(zmat0(ns,:,:),[nm nn])');
%                 xlim([knotsm(1) knotsm(end)]); ylim([knotsn(1) knotsn(end)]);
%                 drawnow;
% 
%             end
% 
%         end
%         
%     end
% 
% %    wvec1 = calcwdn(vmat1(:,1),xmat1(:,1),Gs,Ps,N,bet,sig,kap,rstar)
%     vvec0 = vmat0(:,indm,indn);
%     pvec0 = pmat0(:,indm,indn);
%     xvec0 = xmat0(:,indm,indn);
%     wvec1 = calcwdn(vvec0,pvec0,xvec0,Gs,Ps,N,bet,kap,invlam,sig,rstar);
% 
%     diffout = max(abs(wvec1-wvec0))
%     iterout = iterout + 1
% %    pause
%     wvec0 = wvec1;
%         
% end
% 
% disp(' ');
% toc;
% mmats = mmat0;
% nmats = nmat0;
% 
% eval(sprintf('save wsmat_invlam%6.5fph%1.5fpl%1.2fkap%1.5fN%d.mat vmat0 mmat0 nmat0 xmat0 pmat0 rmat0 zmat0 idmat Gs Ps knotsm knotsn ns nm nn rstar',invlam,ph,pl,kap,N));
% 
% % calculate bounds
% [msvec nsvec mcvec ncvec] = irfs_bound(pcf,pl,knotsm,knotsn,mmats,nmats,mmatc,nmatc);
% 
% disp(' ');
% disp(' the bounds on LMs:'); 
% disp(sprintf(' m: [%10.5f, %10.5f], n:[%10.5f, %10.5f]',mmax,mmin,nmax,nmin));
% disp(' ');
% disp(' the upper and lower values of LMs in IRF:')
% disp(sprintf(' m: (%10.5f, %10.5f), n:(%10.5f, %10.5f)',max(msvec),min(msvec),max(nsvec),min(nsvec)));
% 
% mmax = max(msvec);
% mmin = min(msvec);
% nmax = max(nsvec);
% nmin = min(nsvec);