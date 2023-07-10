clearvars
clc
close all

addpath ./fun 

% Load simulation data, obtained from CREATE-NL
% Kicks start at 1s, and are filtered by the PS model 
load('./data/NL-sim-data.mat')

Ts = diff(Sim.time(1:2)); % Sampling time
Tpulse = 1; % Kicks are given with a T = 1s period

R = 0.0120; % True value of VS3 resistance

% Find first kick
idxStart = find(Sim.time>=1.15,1,'first'); % First time instant after t = 1s
Nper = round(Tpulse/Ts)-1; % Number of samples between two pulses

% Pre-initialize vectors in the for cycle
nn = floor((length(Sim.time)-idxStart)/Nper);
tMod      = zeros(1,nn);
sysZrec   = cell(1,nn);
sysZrecC  = cell(1,nn);
sysZrec0  = cell(1,nn);
sysZrecC0 = cell(1,nn);

k = 1;
for i = idxStart : Nper+1:length(Sim.time)
  
  idx = i : i+500; % Nper = 2000, so this is about 330 samples
  tt = Sim.time(idx);

  tMod(k) = tt(1);
  u = Sim.VS3(idx)-Sim.VS3(idx(1)); % input variation
  y = [Sim.Ivs3(idx)- Sim.Ivs3(idx(1));
       Sim.Zpl(idx) - Sim.Zpl(idx(1))]; % output variation
  
  % Build matrices for discrete-time identification procedure
  T = length(idx)-1;
  X1 = packsignal2(1,T,y);
  X0 = packsignal2(0,T,y);
  U  = packsignal1(0,1,T,u);
  
  % Identify A and B matrices
  BA = X1*pinv([U;X0]);
  B  = BA(:,1);
  A  = BA(:,2:3);
  
  % Fix dcgain: first value (current channel) is determined by R
  G  = (eye(2)-A)\B;
  G0 = [1/R; G(2)]; % Leave second gain unaffected
  
  % Add fictitious regime samples
  NP  = 1;
  W   = 400;
  X00 = [X0 W*repmat(G0,1,NP)];
  X10 = [X1 W*repmat(G0,1,NP)];
  U0  = [U  W*ones(1,NP)];
  
  % Repeat identification
  BA0 = X10*pinv([U0;X00]);
  B0  = BA0(:,1);
  A0  = BA0(:,2:3);
  
  % Store identified system in both discrete-time and continuous-time
  sysZrec{k}  = ss(A,B,eye(2),zeros(2,1),Ts);
  sysZrecC{k} = d2c(sysZrec{k},'zoh'); % Convert from discrete to continuous
  
  sysZrec0{k}  = ss(A0,B0,eye(2),zeros(2,1),Ts);
  sysZrecC0{k} = d2c(sysZrec0{k},'zoh'); % Convert from discrete to continuous
  k = k+1;
end

save .\data\models-gains-constrained sysZrecC sysZrecC0

