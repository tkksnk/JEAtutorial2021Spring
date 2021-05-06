clear all;

m.rstar = 0.75; % pH=0�̂Ƃ��́A����Ԃł̖��ڋ����̒l
m.bet = 1/(1+m.rstar/100); % ������(�I�C���[�������̒���Ԃ��)
m.phi = 5.0;  % �e�C���[�W��(��: ��������iL=0�ɂȂ�Ȃ�)
m.pL = 0.75;  % ��@�̌p���m��
m.sH = m.rstar; % ���H�ł̎��R���q���̒l

% �J���u���[�V����
% yL��piL�̃^�[�Q�b�g��pH=0�̂Ƃ��̃��f���̒l�����킹��悤�ɁAsL��kap�̒l���Z�b�g

m.pH = 0.0; % ��@���N����m��
x0 = [-2.0, 0.01]; % sL��kap�̏����l

% yL��piL�̃^�[�Q�b�g
yLtar = -7.0;
piLtar = -1.0/4;

% �ŏ����֐�(Matlab�̏ꍇfminsearch)��p����
x = fminsearch(@dist,x0,[],m.sH,m.pH,m.pL,m.bet,m.phi,m.rstar,yLtar,piLtar);

% �J���u���[�g�����p�����[�^���Z�b�g
m.sL = x(1); % ���L�ł̎��R���q���̒l
m.kap = x(2); % �t�B���b�v�X�Ȑ��̌X��

m.maxiter = 2000; % �J��Ԃ��񐔂̍ő�l
m.tol = 1e-5; % ���e�덷

%%
tic;
[yvec0 pvec0 rvec0] = ti(m);
toc;

m.pH = 0.025
tic;
[yvec1 pvec1 rvec1] = ti(m);
toc;

%%
xvec = [0 1];

figure;
subplot(231);
plot(xvec,rvec0*4,'k*-','LineWidth',3.0);
hold on;
plot(xvec,rvec1*4,'k*--','LineWidth',3.0);
plot(xvec,[0 0],'r-');
title('�������');
xticks([0 1]);
xticklabels({'H','L'});
set(gca,'Fontsize',12);

subplot(232);
plot(xvec,yvec0,'k*-','LineWidth',3.0);
hold on;
plot(xvec,yvec1,'k*--','LineWidth',3.0);
plot(xvec,[0 0],'r-');
title('�Y�o�M���b�v');
xticks([0 1]);
xticklabels({'H','L'});
set(gca,'Fontsize',12);

subplot(233);
plot(xvec,pvec0*4,'k*-','LineWidth',3.0);
hold on;
plot(xvec,pvec1*4,'k*--','LineWidth',3.0);
plot(xvec,[0 0],'r-');
title('�C���t����');
xticks([0 1]);
xticklabels({'H','L'});
set(gca,'Fontsize',12);
m = legend('p_H=0','p_H=0.025','Location','SouthWest');
m.FontSize = 8;
% saveas(gcf,'simplepf.eps','epsc2');