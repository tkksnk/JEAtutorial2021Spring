% clear all;
% ph = 0.001;
% pl = 0.5;
eval(sprintf('load wsmat_invlam%6.5fph%1.5fpl%1.2fN%d.mat vmat0 mmat0 nmat0 xmat0 pmat0 rmat0 zmat0 idmat Gs Ps knotsm knotsn ns nm nn rstar',invlam,ph,pl,N));
mmin = knotsm(1);
mmax = knotsm(end);
nmin = knotsn(1);
nmax = knotsn(end);

figure
subplot(221)
mesh(knotsm,knotsn,reshape(mmat0(1,:,:),[nm nn])')
title('s=s_h');
zlabel('\phi_{1}');
xlim([mmin mmax]); ylim([nmin nmax]);
% zlim([mmin mmax]);
subplot(222)
mesh(knotsm,knotsn,reshape(mmat0(2,:,:),[nm nn])')
title('s=s_l');
xlim([mmin mmax]); ylim([nmin nmax]);
% zlim([mmin mmax]);
subplot(223)
mesh(knotsm,knotsn,reshape(nmat0(1,:,:),[nm nn])')
xlabel('\phi_{1,-1}'); ylabel('\phi_{2,-1}');
zlabel('\phi_{2}');
xlim([mmin mmax]); ylim([nmin nmax]);
% zlim([nmin nmax]);
subplot(224)
mesh(knotsm,knotsn,reshape(nmat0(2,:,:),[nm nn])')
xlabel('\phi_{1,-1}'); ylabel('\phi_{2,-1}');
xlim([mmin mmax]); ylim([nmin nmax]);
% zlim([nmin nmax]);

%eval(sprintf('print -depsc2 pfs1_ph%1.5fpl%1.2f.eps',ph,pl));

figure
subplot(321)
mesh(knotsm,knotsn,reshape(xmat0(1,:,:),[nm nn])')
title('s=s_h');
zlabel('x');
xlim([mmin mmax]); ylim([nmin nmax]);
subplot(322)
mesh(knotsm,knotsn,reshape(xmat0(2,:,:),[nm nn])')
title('s=s_l');
xlim([mmin mmax]); ylim([nmin nmax]);
subplot(323)
mesh(knotsm,knotsn,reshape(pmat0(1,:,:),[nm nn])')
zlabel('\pi');
xlim([mmin mmax]); ylim([nmin nmax]);
subplot(324)
mesh(knotsm,knotsn,reshape(pmat0(2,:,:),[nm nn])')
xlim([mmin mmax]); ylim([nmin nmax]);
subplot(325)
mesh(knotsm,knotsn,reshape(rmat0(1,:,:),[nm nn])')
xlabel('\phi_{1,-1}'); ylabel('\phi_{2,-1}');
zlabel('r');
xlim([mmin mmax]); ylim([nmin nmax]);
subplot(326)
mesh(knotsm,knotsn,reshape(rmat0(2,:,:),[nm nn])')
xlabel('\phi_{1,-1}'); ylabel('\phi_{2,-1}');
xlim([mmin mmax]); ylim([nmin nmax]);

%eval(sprintf('print -depsc2 pfs2_ph%1.5fpl%1.2f.eps',ph,pl));

figure
subplot(321)
mesh(knotsm,knotsn,reshape(zmat0(1,:,:),[nm nn])')
title('s=s_h');
zlabel('z');
xlim([mmin mmax]); ylim([nmin nmax]);
zlim([0 1]);
subplot(322)
mesh(knotsm,knotsn,reshape(zmat0(2,:,:),[nm nn])')
title('s=s_l');
xlim([mmin mmax]); ylim([nmin nmax]);
zlim([0 1]);
subplot(323)
mesh(knotsm,knotsn,reshape(vmat0(1,:,:),[nm nn])')
xlabel('\phi_{1,-1}'); ylabel('\phi_{2,-1}');
zlabel('v');
xlim([mmin mmax]); ylim([nmin nmax]);
subplot(324)
mesh(knotsm,knotsn,reshape(vmat0(2,:,:),[nm nn])')
xlabel('\phi_{1,-1}'); ylabel('\phi_{2,-1}');
xlim([mmin mmax]); ylim([nmin nmax]);
subplot(325)
mesh(knotsm,knotsn,reshape(idmat(1,:,:),[nm nn])')
xlabel('\phi_{1,-1}'); ylabel('\phi_{2,-1}');
zlabel('I'); 
xlim([mmin mmax]); ylim([nmin nmax]);
zlim([1 4]);
subplot(326)
mesh(knotsm,knotsn,reshape(idmat(2,:,:),[nm nn])')
xlabel('\phi_{1,-1}'); ylabel('\phi_{2,-1}');
xlim([mmin mmax]); ylim([nmin nmax]);
zlim([1 4]);

%eval(sprintf('print -depsc2 pfs3_ph%1.5fpl%1.2f.eps',ph,pl));
