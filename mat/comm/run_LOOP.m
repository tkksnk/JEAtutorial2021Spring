% run_LOOP.m
% adaptive grid points based on IRFs.
clear all;
% 
pcf = 0.001;  % cutoff probability for calculating bounds: 
% the likelihood of the low state duration T1 is pl^T1.
% taking as given the cutoff probability pcf, i.e., pl^T1>pcf, T1 is 
% T1 = ceil(log(pcf)/log(pl))
xdtar = -7.0; % crisis output under the discretionary policy
pdtar = -1.0/4;

invlam = 16;  % > 105 does not work
ph = 0.005;   % crisis probability: ph = 1e-8 yields MPE as sustainable policy
pl = 0.5;     % crisis duration: > 0.6 may not work 
N = 80;      % length of punishment periods

% bounds on Lagrange multipliers: first we run the model with a wider
% range, and then reset the range based on simulation and re-run the model.
mmax = 0.0;
mmin = -200.0;
nmax = 50.0;
nmin = 0.0;

[mmax mmin nmax nmin sl kap] = mainsp_LOOP(xdtar,pdtar,ph,pl,invlam,mmin,mmax,nmin,nmax,pcf,N);

[mmax mmin nmax nmin sl kap] = mainsp_LOOP(xdtar,pdtar,ph,pl,invlam,mmin,mmax,nmin,nmax,pcf,N);

% plot irfs
plotirfs;
eval(sprintf('save irfmat_ph%1.5fpl%1.2fN%d.mat vsvec nsvec xsvec rsvec zsvec vcvec ncvec xcvec rcvec vdvec xdvec rdvec',ph,pl,N));
%plotirfs2_092717