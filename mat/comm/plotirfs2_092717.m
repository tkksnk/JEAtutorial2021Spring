% plotirf.m: plot impulse responce functions

eval(sprintf('load wsmat_invlam%6.5fph%1.5fpl%1.2fN%d.mat vmat0 mmat0 nmat0 xmat0 pmat0 rmat0 zmat0 idmat Gs Ps knotsm knotsn ns nm nn rstar',invlam,ph,pl,N));
vmats = vmat0;
mmats = mmat0;
nmats = nmat0;
xmats = xmat0;
pmats = pmat0;
rmats = rmat0;
zmats = zmat0;
eval(sprintf('load wcmat_invlam%6.5fph%1.5fpl%1.2fN%d.mat vmat0 mmat0 nmat0 xmat0 pmat0 rmat0 Gs Ps knotsm knotsn ns nm nn rstar',invlam,ph,pl,N));
vmatc = vmat0;
mmatc = mmat0;
nmatc = nmat0;
xmatc = xmat0;
pmatc = pmat0;
rmatc = rmat0;
eval(sprintf('load wdmat_invlam%6.5fph%1.5fpl%1.2fN%d.mat vvec0 xvec0 pvec0 rvec0 Gs Ps ns rstar',invlam,ph,pl,N));
vvecd = vvec0;
xvecd = xvec0;
pvecd = pvec0;
rvecd = rvec0;


T1 = ceil(log(pcf)/log(pl));
T = 2*T1;
%T = 20;
%T1 = 8;

% initial values
vsvec = zeros(T+2,1);
msvec = zeros(T+2,1);
nsvec = zeros(T+2,1);
xsvec = zeros(T+2,1);
psvec = zeros(T+2,1);
rsvec = zeros(T+2,1);
zsvec = ones(T+2,1);
vcvec = zeros(T+2,1);
mcvec = zeros(T+2,1);
ncvec = zeros(T+2,1);
xcvec = zeros(T+2,1);
pcvec = zeros(T+2,1);
rcvec = zeros(T+2,1);
vdvec = zeros(T+2,1);
xdvec = zeros(T+2,1);
pdvec = zeros(T+2,1);
rdvec = zeros(T+2,1);

% compute initial value (stochastic steady state)
% ivec = ones(1000,1);
% 
% for t=1:1000
% 
%     msvec(t+1) = intf2(knotsm,knotsn,reshape(mmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
%     nsvec(t+1) = intf2(knotsm,knotsn,reshape(nmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
%     vsvec(t+1) = intf2(knotsm,knotsn,reshape(vmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
%     xsvec(t+1) = intf2(knotsm,knotsn,reshape(xmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
%     psvec(t+1) = intf2(knotsm,knotsn,reshape(pmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
%     rsvec(t+1) = intf2(knotsm,knotsn,reshape(rmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
% 
%     mcvec(t+1) = intf2(knotsm,knotsn,reshape(mmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
%     ncvec(t+1) = intf2(knotsm,knotsn,reshape(nmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
%     vcvec(t+1) = intf2(knotsm,knotsn,reshape(vmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
%     xcvec(t+1) = intf2(knotsm,knotsn,reshape(xmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
%     pcvec(t+1) = intf2(knotsm,knotsn,reshape(pmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
%     rcvec(t+1) = intf2(knotsm,knotsn,reshape(rmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
%    
% end
% 
% msvec(1) = msvec(end);
% nsvec(1) = nsvec(end);
% vsvec(1) = vsvec(end);
% xsvec(1) = xsvec(end);
% psvec(1) = psvec(end);
% rsvec(1) = rsvec(end);
% 
% mcvec(1) = mcvec(end);
% ncvec(1) = ncvec(end);
% vcvec(1) = vcvec(end);
% xcvec(1) = xcvec(end);
% pcvec(1) = pcvec(end);
% rcvec(1) = rcvec(end);
% 
% vdvec(1) = vvecd(1);
% xdvec(1) = xvecd(1);
% pdvec(1) = pvecd(1);
% rdvec(1) = rvecd(1);
esmat = zeros(T+2,3);
ecmat = zeros(T+2,3);

ivec = ones(T+1,1);
ivec(1:T1) = 2;

for t=1:T+1

    msvec(t+1) = intf2(knotsm,knotsn,reshape(mmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
    nsvec(t+1) = intf2(knotsm,knotsn,reshape(nmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
    vsvec(t+1) = intf2(knotsm,knotsn,reshape(vmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
    xsvec(t+1) = intf2(knotsm,knotsn,reshape(xmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
    psvec(t+1) = intf2(knotsm,knotsn,reshape(pmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
    rsvec(t+1) = intf2(knotsm,knotsn,reshape(rmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
    zsvec(t+1) = intf2(knotsm,knotsn,reshape(zmats(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
    idvec(t+1) = intf2(knotsm,knotsn,reshape(idmat(ivec(t),:,:),[nm nn]),msvec(t),nsvec(t));
        
    mcvec(t+1) = intf2(knotsm,knotsn,reshape(mmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
    ncvec(t+1) = intf2(knotsm,knotsn,reshape(nmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
    vcvec(t+1) = intf2(knotsm,knotsn,reshape(vmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
    xcvec(t+1) = intf2(knotsm,knotsn,reshape(xmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
    pcvec(t+1) = intf2(knotsm,knotsn,reshape(pmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
    rcvec(t+1) = intf2(knotsm,knotsn,reshape(rmatc(ivec(t),:,:),[nm nn]),mcvec(t),ncvec(t));
    
    vdvec(t+1) = vvecd(ivec(t));
    xdvec(t+1) = xvecd(ivec(t));
    pdvec(t+1) = pvecd(ivec(t));
    rdvec(t+1) = rvecd(ivec(t));
    
%     % short-run gain and long-run loss
    usvec(t+1) = -(xsvec(t+1)^2 + invlam*psvec(t+1)^2);
    ucvec(t+1) = -(xcvec(t+1)^2 + invlam*pcvec(t+1)^2);
    udvec(t+1) = -(xdvec(t+1)^2 + invlam*pdvec(t+1)^2);
%     
%     evsvec(t+1) = vsvec(t+1) - usvec(t+1);
%     evcvec(t+1) = vcvec(t+1) - ucvec(t+1);
%     evdvec(t+1) = vdvec(t+1) - udvec(t+1);
%     
%     SRgsvec(t+1) = udvec(t+1) - usvec(t+1);
%     SRgcvec(t+1) = udvec(t+1) - ucvec(t+1);
%     LRgsvec(t+1) = evsvec(t+1) - evdvec(t+1);
%     LRgcvec(t+1) = evcvec(t+1) - evdvec(t+1);    
%    
%     % Euler equation error. 11/22/16: need to be fixed?
%     vcond = Ps(ivec(t),1)*vmats(1,:,:) + Ps(ivec(t),2)*vmats(2,:,:); 
%     xcond = Ps(ivec(t),1)*xmats(1,:,:) + Ps(ivec(t),2)*xmats(2,:,:); 
%     pcond = Ps(ivec(t),1)*pmats(1,:,:) + Ps(ivec(t),2)*pmats(2,:,:); 
%     vcond = reshape(vcond,[nm nn]);
%     pcond = reshape(pcond,[nm nn]);
%     xcond = reshape(xcond,[nm nn]);
%     esvec = eulererr(Gs(ivec(t)),msvec(t+1),nsvec(t+1),vsvec(t+1),xsvec(t+1),psvec(t+1),rsvec(t+1),...
%         knotsm,knotsn,vcond,xcond,pcond,bet,kap,invlam,sig);
%     esmat(t+1,1) = esvec(1);
%     esmat(t+1,2) = esvec(2);
%     esmat(t+1,3) = esvec(3);
% 
%     vcond = Ps(ivec(t),1)*vmatc(1,:,:) + Ps(ivec(t),2)*vmatc(2,:,:); 
%     xcond = Ps(ivec(t),1)*xmatc(1,:,:) + Ps(ivec(t),2)*xmatc(2,:,:); 
%     pcond = Ps(ivec(t),1)*pmatc(1,:,:) + Ps(ivec(t),2)*pmatc(2,:,:); 
%     vcond = reshape(vcond,[nm nn]);
%     pcond = reshape(pcond,[nm nn]);
%     xcond = reshape(xcond,[nm nn]);
%     ecvec = eulererr(Gs(ivec(t)),mcvec(t+1),ncvec(t+1),vcvec(t+1),xcvec(t+1),pcvec(t+1),rcvec(t+1),...
%         knotsm,knotsn,vcond,xcond,pcond,bet,kap,invlam,sig);
%     ecmat(t+1,1) = ecvec(1);
%     ecmat(t+1,2) = ecvec(2);
%     ecmat(t+1,3) = ecvec(3);

end


% load parms rstar;
rstar = 0.0075*100;
rsvec = rsvec + rstar;
rcvec = rcvec + rstar;
rdvec = rdvec + rstar;

% plot figures
figure;

subplot(221);
plot([-1:T],xcvec(1:T+2),'gx-','LineWidth',2.0);
hold on;
plot([-1:T],xdvec,'r+-','LineWidth',2.0);
plot([-1:T],xsvec(1:T+2),'bo-','LineWidth',2.0);
title('Output gap (%)');
xlim([-1 T]);

subplot(222);
plot([-1:T],rcvec(1:T+2)*4,'gx-','LineWidth',2.0);
hold on;
plot([-1:T],rdvec*4,'r+-','LineWidth',2.0);
plot([-1:T],rsvec(1:T+2)*4,'bo-','LineWidth',2.0);
title('Nom. interest rate (%)');
xlim([-1 T]);
% legend('Sustainable','Commitment','Discretion','Location','SouthEast');

subplot(223);
plot([-1:T],vcvec(1:T+2),'gx-','LineWidth',2.0);
hold on;
plot([-1:T],vdvec,'r+-','LineWidth',2.0);
plot([-1:T],vsvec(1:T+2),'bo-','LineWidth',2.0);
title('Value');
xlim([-1 T]);

subplot(224);
plot([-1:T],pcvec(1:T+2)*4,'gx-','LineWidth',2.0);
hold on;
plot([-1:T],pdvec*4,'r+-','LineWidth',2.0);
plot([-1:T],psvec(1:T+2)*4,'bo-','LineWidth',2.0);
title('Inflation rate (%)');
xlim([-1 T]);

save irfs2_ph0.002500.mat xcvec xdvec xsvec rcvec rdvec rsvec vcvec vdvec vsvec pcvec pdvec psvec;

% subplot(425);
% plot([-1:T],ucvec(1:T+2),'gx-','LineWidth',2.0);
% hold on;
% plot([-1:T],udvec(1:T+2),'r+-','LineWidth',2.0);
% plot([-1:T],usvec(1:T+2),'bo-','LineWidth',2.0);
% title('Utility');
% xlim([-1 T]);
% 
% subplot(426);
% plot([-1:T],zsvec(1:T+2),'gx-','LineWidth',2.0);
% hold on;
% plot([-1:T],zsvec(1:T+2),'bo-','LineWidth',2.0);
% title('z');
% xlim([-1 T]);
% 
% subplot(427);
% plot([-1:T],mcvec(1:T+2),'gx-','LineWidth',2.0);
% hold on;
% plot([-1:T],msvec(1:T+2),'bo-','LineWidth',2.0);
% %title('\phi_{1t}');
% title('LM on PC');
% xlim([-1 T]);
% 
% subplot(428);
% plot([-1:T],ncvec(1:T+2),'gx-','LineWidth',2.0);
% hold on;
% plot([-1:T],nsvec(1:T+2),'bo-','LineWidth',2.0);
% %title('\phi_{2t}');
% title('LM on Euler');
% xlim([-1 T]);

%eval(sprintf('print -depsc2 irf%d_ph%1.5fpl%1.2f.eps',T1,ph,pl));

% figure;
% subplot(211);
% plot([-1:T],vcvec(1:T+2)-vdvec(1:T+2),'g-','LineWidth',2.0);
% hold on;
% plot([-1:T],vsvec(1:T+2)-vdvec(1:T+2),'b-','LineWidth',2.0);
% plot([-1 T],[0 0],'r-');
% xlim([-1 T]);
% 
% subplot(212);
% plot([-1:T],SRgcvec(1:T+2),'g-','LineWidth',2.0);
% hold on;
% plot([-1:T],SRgsvec(1:T+2),'b-','LineWidth',2.0);
% plot([-1:T],LRgcvec(1:T+2),'g--','LineWidth',2.0);
% plot([-1:T],LRgsvec(1:T+2),'b--','LineWidth',2.0);
% % legend('Short-run gain Ramsey','Short-run gain Sustainable', ...
% % 'Long-run loss Ramsey','Long-run loss Sustainable');
% plot([-1 T],[0 0],'r-');
% xlim([-1 T]);
