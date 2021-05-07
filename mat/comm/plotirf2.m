clear all;

load irfs2_ph0.00050.mat xcvec xdvec xsvec rcvec rdvec rsvec vcvec vdvec vsvec pcvec pdvec psvec;
xcvec1 = xcvec;
xdvec1 = xdvec;
xsvec1 = xsvec;
rcvec1 = rcvec;
rdvec1 = rdvec;
rsvec1 = rsvec;
vcvec1 = vcvec;
vdvec1 = vdvec;
vsvec1 = vsvec;
pcvec1 = pcvec;
pdvec1 = pdvec;
psvec1 = psvec;
load irfs2_ph0.00100.mat xcvec xdvec xsvec rcvec rdvec rsvec vcvec vdvec vsvec pcvec pdvec psvec;
xcvec2 = xcvec;
xdvec2 = xdvec;
xsvec2 = xsvec;
rcvec2 = rcvec;
rdvec2 = rdvec;
rsvec2 = rsvec;
vcvec2 = vcvec;
vdvec2 = vdvec;
vsvec2 = vsvec;
pcvec2 = pcvec;
pdvec2 = pdvec;
psvec2 = psvec;
load irfs2_ph0.002500.mat xcvec xdvec xsvec rcvec rdvec rsvec vcvec vdvec vsvec pcvec pdvec psvec;
xcvec3 = xcvec;
xdvec3 = xdvec;
xsvec3 = xsvec;
rcvec3 = rcvec;
rdvec3 = rdvec;
rsvec3 = rsvec;
vcvec3 = vcvec;
vdvec3 = vdvec;
vsvec3 = vsvec;
pcvec3 = pcvec;
pdvec3 = pdvec;
psvec3 = psvec;

% plot figures
T = 16;
figure;

subplot(221);
plot([-1:T],xcvec1(1:T+2),'g-','LineWidth',2.0);
hold on;
plot([-1:T],xsvec3(1:T+2),'b:','LineWidth',2.0);
plot([-1:T],xsvec2(1:T+2),'b--','LineWidth',2.0);
plot([-1:T],xsvec1(1:T+2),'b-','LineWidth',2.0);
plot([-1:T],xdvec1,'r-','LineWidth',2.0);
plot([-1:T],xcvec2(1:T+2),'g--','LineWidth',2.0);
%hold on;
plot([-1:T],xdvec2,'r--','LineWidth',2.0);
plot([-1:T],xcvec3(1:T+2),'g:','LineWidth',2.0);
%hold on;
plot([-1:T],xdvec3,'r:','LineWidth',2.0);
title('Output gap (%)');
xlim([-1 T]);
yticks([-5,-2.5,0,2.5,5]);
legend({'Ramsey','OSP:100p_{h}=0.25','OSP:100p_{h}=0.1','OSP:100p_{h}=0.05','Markov-Perfect'},'FontSize',10,'Location','NorthWest')

subplot(223);
plot([-1:T],rcvec1(1:T+2)*4,'g-','LineWidth',2.0);
hold on;
plot([-1:T],rdvec1*4,'r-','LineWidth',2.0);
plot([-1:T],rsvec1(1:T+2)*4,'b-','LineWidth',2.0);
plot([-1:T],rcvec2(1:T+2)*4,'g--','LineWidth',2.0);
%hold on;
plot([-1:T],rdvec2*4,'r--','LineWidth',2.0);
plot([-1:T],rsvec2(1:T+2)*4,'b--','LineWidth',2.0);
plot([-1:T],rcvec3(1:T+2)*4,'g:','LineWidth',2.0);
%hold on;
plot([-1:T],rdvec3*4,'r:','LineWidth',2.0);
plot([-1:T],rsvec3(1:T+2)*4,'b:','LineWidth',2.0);
title('Policy rate (%)');
xlim([-1 T]);
ylim([0 3.5]);
yticks([0,1,2,3]);
% legend('Sustainable','Commitment','Discretion','Location','SouthEast');

% subplot(224);
% plot([-1:T],vsvec1(1:T+2),'b-','LineWidth',2.0);
% hold on;
% plot([-1:T],vsvec2(1:T+2),'b--','LineWidth',2.0);
% plot([-1:T],vsvec3(1:T+2),'b:','LineWidth',2.0);
% plot([-1:T],vcvec1(1:T+2),'g-','LineWidth',2.0);
% plot([-1:T],vdvec1,'r-','LineWidth',2.0);
% plot([-1:T],vcvec2(1:T+2),'g--','LineWidth',2.0);
% %hold on;
% plot([-1:T],vdvec2,'r--','LineWidth',2.0);
% plot([-1:T],vcvec3(1:T+2),'g:','LineWidth',2.0);
% %hold on;
% plot([-1:T],vdvec3,'r:','LineWidth',2.0);
% title('Value');
% xlim([-1 T]);

subplot(222);
plot([-1:T],pcvec1(1:T+2)*4,'g-','LineWidth',2.0);
hold on;
plot([-1:T],psvec1(1:T+2)*4,'b-','LineWidth',2.0);
plot([-1:T],psvec2(1:T+2)*4,'b--','LineWidth',2.0);
plot([-1:T],psvec3(1:T+2)*4,'b:','LineWidth',2.0);
plot([-1:T],pdvec1*4,'r-','LineWidth',2.0);
plot([-1:T],pcvec2(1:T+2)*4,'g--','LineWidth',2.0);
%hold on;
plot([-1:T],pdvec2*4,'r--','LineWidth',2.0);
plot([-1:T],pcvec3(1:T+2)*4,'g:','LineWidth',2.0);
%hold on;
plot([-1:T],pdvec3*4,'r:','LineWidth',2.0);
title('Inflation rate (%)');
xlim([-1 T]);

set(gca,'FontWeight','normal');
