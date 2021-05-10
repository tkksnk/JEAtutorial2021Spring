clear all;
load pf.mat;

N = size(ymat0,1);
Ns = size(ymat0,2);
mmin = mgrid(1);
mmax = mgrid(end);

% simulation
T = 40;
T1 = 8;
svec = ones(T+1,1);
yvec = zeros(T+1,1);
pvec = zeros(T+1,1);
ivec = zeros(T+1,1);
mvec = zeros(T+1,1);
nvec = zeros(T+1,1);

svec(1:T1) = 2;

options = optimset('TolFun',1e-10);

for t=1:T

    is = svec(t);
    s0 = Gs(is);
    mpast = mvec(t);

    eyvec = zeros(N,1);
    epvec = zeros(N,1);
    for js = 1:Ns

        eyvec = eyvec + Ps(is,js)*ymat0(:,js);
        epvec = epvec + Ps(is,js)*pmat0(:,js);

    end
    
    % solve for m0
    x = fzero(@rw_eqm,mpast,options,mpast,mgrid,epvec,rstar,phi,zeta);
    %x = golden('rw_eqm',mmin,mmax,1e-10,mpast,mgrid,epvec,rstar,phi,zeta);
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
    
    yvec(t+1) = y0;
    pvec(t+1) = p0;
    ivec(t+1) = i0;
    mvec(t+1) = m0;
    nvec(t+1) = is0+mpast;

end

yvec = yvec(2:T+1);
pvec = pvec(2:T+1);
ivec = ivec(2:T+1);
mvec = mvec(2:T+1);
nvec = nvec(2:T+1);

if (zeta>0)
    temp = abs(mgrid-mmat0(:,2));
    id = find(min(temp)==temp);
else
    id = 1;
end
mgrid(id)

figure;
subplot(121);
plot(mgrid,mmat0(:,1),'b-','LineWidth',2.0);
hold on;
plot(mgrid,mgrid,'k:');
plot(mvec(T1:end-1),mvec(T1+1:end),'m*');
xlim([mgrid(1) mgrid(end)]);
ylim([mgrid(1) mgrid(end)]);
title('s_{H}','FontWeight','Normal');
xlabel('m_{t-1}');
ylabel('m_{t}','FontWeight','Normal')
subplot(122);
plot(mgrid,mmat0(:,2),'b-','LineWidth',2.0);
hold on;
plot(mgrid,mgrid,'k:');
plot([0;mvec(1:T1-1)],mvec(1:T1),'m*');
xlim([mgrid(1) mgrid(end)]);
ylim_ = ylim;
plot([mgrid(id) mgrid(id)],[ylim_(1) ylim_(2)],'r-');
title('s_{L}','FontWeight','Normal');
xlabel('m_{t-1}');
ylabel('m_{t}','FontWeight','Normal')

T = 20;
figure;
% subplot(231);
subplot(232);
plot([1:T],mvec(1:T),'k-','LineWidth',2.0);
title('m_{t}','FontWeight','Normal');
xlim([1 T]);
xticks([1 10 20]);
subplot(231);
plot([1:T],Gs(svec(1:T)),'r-','LineWidth',2.0);
title('s_{t}','FontWeight','Normal');
xlim([1 T]);
xticks([1 10 20]);
% subplot(232);
% plot([1:T],nvec(1:T)*4,'k-','LineWidth',2.0);
% title('r_{n,t}^{*}','FontWeight','Normal');
% xlim([1 T]);
% xticks([1 10 20]);
% subplot(233);
% plot([1:T],ivec(1:T)*4,'k-','LineWidth',2.0);
% title('r_{n,t}','FontWeight','Normal');
% xlim([1 T]);
% xticks([1 10 20]);

% figure;
% subplot(231);
% plot(mgrid,mmat0(:,1),'b-','LineWidth',2.0);
% hold on;
% plot(mgrid,mgrid,'k:');
% plot(mvec(T1:end-1),mvec(T1+1:end),'m*');
% xlim([mgrid(1) mgrid(end)]);
% ylim([mgrid(1) mgrid(end)]);
% ylabel('s_{H}');
% title('FG term, m_{t}','FontWeight','Normal')
% subplot(234);
% plot(mgrid,mmat0(:,2),'b-','LineWidth',2.0);
% hold on;
% plot(mgrid,mgrid,'k:');
% plot([0;mvec(1:T1-1)],mvec(1:T1),'m*');
% xlim([mgrid(1) mgrid(end)]);
% ylim_ = ylim;
% plot([mgrid(id) mgrid(id)],[ylim_(1) ylim_(2)],'r-');
% xlabel('m_{t-1}');
% ylabel('s_{L}');
% 
% subplot(232);
% plot(mgrid,nmat0(:,1),'b-','LineWidth',2.0);
% xlim([mgrid(1) mgrid(end)]);
% ylim([-4 1]);
% ylim_ = ylim;
% hold on;
% plot(mvec(T1:end-1),nvec(T1+1:end),'m*');
% plot([mgrid(id) mgrid(id)],[ylim_(1) ylim_(2)],'r-');
% plot([mgrid(1) mgrid(end)],[0,0],'k-');
% title('Nortional rate, R^{*}_{t}-m_{t-1}','FontWeight','Normal')
% subplot(235);
% plot(mgrid,nmat0(:,2),'b-','LineWidth',2.0);
% xlim([mgrid(1) mgrid(end)]);
% ylim([-4 1]);
% ylim_ = ylim;
% hold on;
% plot([0;mvec(1:T1-1)],nvec(1:T1),'m*');
% plot([mgrid(id) mgrid(id)],[ylim_(1) ylim_(2)],'r-');
% plot([mgrid(1) mgrid(end)],[0,0],'k-');
% xlabel('m_{t-1}')
% 
% subplot(233);
% plot(mgrid,imat0(:,1),'b-','LineWidth',2.0);
% xlim([mgrid(1) mgrid(end)]);
% ylim([0 1]);
% ylim_ = ylim;
% hold on;
% plot(mvec(T1:end-1),ivec(T1+1:end),'m*');
% plot([mgrid(id) mgrid(id)],[ylim_(1) ylim_(2)],'r-');
% title('Nominal rate, R_{t}','FontWeight','Normal')
% subplot(236);
% plot(mgrid,imat0(:,2),'b-','LineWidth',2.0);
% xlim([mgrid(1) mgrid(end)]);
% ylim([0 1]);
% ylim_ = ylim;
% hold on;
% plot([0;mvec(1:T1-1)],ivec(1:T1),'m*');
% plot([mgrid(id) mgrid(id)],[ylim_(1) ylim_(2)],'r-');
% % plot([mgrid(1) mgrid(end)],[0,0],'k-');
% xlabel('m_{t-1}')
% 
% Tshow = 20;
% figure;
% subplot(232);
% plot([1:Tshow],yvec(1:Tshow),'o-');
% xlim([1 Tshow]);
% title('Output');
% subplot(233);
% plot([1:Tshow],pvec(1:Tshow)*4,'o-');
% xlim([1 Tshow]);
% title('Inflation');
% subplot(231);
% plot([1:Tshow],ivec(1:Tshow)*4,'o-');
% xlim([1 Tshow]);
% title('Policy rate');
% 
% save rwirf2.5.mat yvec pvec ivec;
% subplot(234);
% plot([1:Tshow],mvec(1:Tshow),'o-');
% xlim([1 Tshow]);
% title('FG term');
% subplot(235);
% plot([1:Tshow],nvec(1:Tshow),'o-');
% hold on;
% plot([1:Tshow],nvec(1:Tshow)+mvec(1:Tshow),'o-');
% plot([1:Tshow],ivec(1:Tshow),'o-');
% xlim([1 Tshow]);
% title('Nortional rate');
% subplot(236);
% plot([1:Tshow],Gs(svec(1:Tshow)),'ro-');
% xlim([1 Tshow]);
% title('Natural rate');

% save irf.mat yvec pvec ivec mvec;
