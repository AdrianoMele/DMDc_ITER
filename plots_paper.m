clearvars
close all
clc

addpath ./fun
%% Load input data
load('./data/Eq_RD_nwController_2d4s_CL.mat') % Plasma linearized model
[A_contr,B_contr,C_contr,D_contr] = initVertPlasmaModel(LinearModel);
sysFull = ss(A_contr,B_contr,C_contr,D_contr);

load('./data/models-gains-constrained.mat')
iModel = 3; 
sysz = minreal(sysZrecC0{iModel});
% sysZrec{k}  = ss(A,B,eye(2),zeros(2,1),Ts);
% sysz = d2c(sysz,'zoh'); % Convert from discrete to continuous

close all
figure('Position',[280 400 700 360])
subplot(121)
bode(minreal(sysFull(1)));
hold on, grid on
bode(minreal(sysz(1)));
title('$I_{VS3}$','Interpreter','latex','fontsize',18)

subplot(122)
bode(minreal(sysFull(2)));
hold on, grid on
bode(minreal(sysz(2)));
title('$z_p$','Interpreter','latex','fontsize',18)

%% Assemble systems for plot
filt  = tf([1 0],[1e-3 1]);
diagn = tf(1, [3e-3, 1]);
psup  = tf(1, [7.5e-3 1]);
psup.InputDelay = 2.5e-3;

psup = pade(psup,2);

% Design controller on nominal system
[K,Reg] = designVS(sysz);

% Add power supplies and derivative on z channel
sysZ = [1 0; 0 filt]*sysz*psup; 

% Simplified open-loop system from 2x2 identified model
sysZ.outputname = {'I_{VS3}','zdot'};
sysZ.inputname = {'V_{VS3}'};
sysOL = series(sysZ,K,'name');

% Full open-loop system from CREATE-L linearized model
sysFullOL = [1 0; 0 filt]*sysFull*psup;
sysFullOL = minreal(sysFullOL);
sysFullOL.inputname  = 'V_{VS3}';
sysFullOL.outputname = {'I_{VS3}','zdot'};
sysFullOL = series(sysFullOL,K,'name');

% With dynamic compensator
sysOL2 = series(Reg,sysOL,'name');
sysFullOL2 = series(Reg,sysFullOL,'name');

%% Model parameters
A = sysz.a;
B = sysz.b;

RLinv = -A(1,1)+A(1,2)*A(2,1)/A(2,2);
LL = 1/(B(1) - A(1,2)*B(2)/A(2,2));
RR = RLinv*LL;

% Display system in zpk form
zpk(tf(sysz))

%% Bode plots
figure
bode(sysOL)
hold on, grid on
bode(sysFullOL,'--r')
legend({'identified','CREATE-L'},'location','best','FontSize',12)

figure
bode(sysOL)
hold on, grid on
bode(sysFullOL,'--r')
bode(sysOL2,':')
bode(sysFullOL2,'.-r')
legend({'identified, static','CREATE-L, static', 'identified, with lag network','CREATE-L, with lag network'},'location','best','FontSize',12)

%% Nichols plot
figure
nicholsplot(sysOL2)
hold on, grid on
nicholsplot(sysFullOL2,'--r')
legend({'identified','CREATE-L'},'location','best','FontSize',12)
