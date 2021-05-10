clear all;

load ../comm/comm.mat;
% load ./rwirf2.5.mat;
% yvec1 = yvec;
% pvec1 = pvec;
% ivec1 = ivec;
load ./rwirf5.mat;
yvec2 = yvec;
pvec2 = pvec;
ivec2 = ivec;
load ./rwirf100.mat;
yvec3 = yvec;
pvec3 = pvec;
ivec3 = ivec;

T = 20;
% plot figures
figure;

subplot(231);
plot([1:T],xcvec(1:T),'k-','LineWidth',2.0);
hold on;
% plot([1:T],yvec1(1:T),'k-.','LineWidth',2.0);
plot([1:T],yvec2(1:T),'k--','LineWidth',2.0);
plot([1:T],yvec3(1:T),'k:','LineWidth',2.0);
plot([1:T],xdvec(1:T),'r-','LineWidth',2.0);
title('Output gap','FontWeight','Normal');
xlim([1 T]);
xticks([1 10 20]);

subplot(233);
plot([1:T],rcvec(1:T)*4,'k-','LineWidth',2.0);
hold on;
% plot([1:T],ivec1(1:T)*4,'k-.','LineWidth',2.0);
plot([1:T],ivec2(1:T)*4,'k--','LineWidth',2.0);
plot([1:T],ivec3(1:T)*4,'k:','LineWidth',2.0);
plot([1:T],rdvec(1:T)*4,'r-','LineWidth',2.0);
title('Policy rate','FontWeight','Normal');
xlim([1 T]);
xticks([1 10 20]);
% legend('Sustainable','Commitment','Discretion','Location','SouthEast');

% subplot(233);
% plot([-1:T],vcvec(1:T+2),'gx-','LineWidth',2.0);
% hold on;
% plot([-1:T],vdvec,'r+-','LineWidth',2.0);
% % plot([-1:T],vsvec(1:T+2),'bo-','LineWidth',2.0);
% title('Value');
% xlim([-1 T]);

subplot(232);
plot([1:T],pcvec(1:T)*4,'k-','LineWidth',2.0);
hold on;
% plot([1:T],pvec1(1:T)*4,'k-.','LineWidth',2.0);
plot([1:T],pvec2(1:T)*4,'k--','LineWidth',2.0);
plot([1:T],pvec3(1:T)*4,'k:','LineWidth',2.0);
plot([1:T],pdvec(1:T)*4,'r-','LineWidth',2.0);
title('Inflation','FontWeight','Normal');
% legend('OCP','RW \phi=2.5','RW \phi=5.0','RW \phi=100.0','ODP');
xlim([1 T]);
xticks([1 10 20]);