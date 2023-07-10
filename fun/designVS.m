function [K,R] = designVS(sys)
 
filt  = tf([1 0],[1e-3 1]);
diagn = tf(1, [3e-3, 1]);
psup  = tf(1, [7.5e-3 1]);
psup.InputDelay = 2.5e-3;

% sysZ = diagn*[1 0; 0 filt]*sys*psup;
sysZ = [1 0; 0 filt]*sys*psup;

sysZ.inputname  = 'V_{VS3}';
sysZ.outputname = {'I_{VS3}','zdot'};

%% Design controller
% Extract zeros and dc gains
[z,~,k] = zpkdata(sys);

Gi = k(1);
Gz = k(2);

gamma1 = z{1}; % this is the zero of the VVS->IVS channel, which is close to the growth rate gamma
z1     = -z{2}; % this is the zero of the VVS->zdot channel, changed in sign

%% Design procedure
% We assume that the system has the following form:
%  GI(s) = Gi (s-gamma1)/(s-gamma)(s+p2)
%  GZ(s) = Gz s(s+z1)/(s-gamma)(s+p2)
%
% Note that, in the transfer function for the IVS channel, the nonminimum
% phase zero gamma1 is close to the growth rate gamma. Physically, this
% depends on the fact that the unstable mode is not very detectable from
% the measurement of the current. On the other side, there is a zero in the
% zdot channel too, which depends on the direct feedthrough of the voltage
% on the zdot channel (D~=0).
% The poles of the two transfer functions are obviously equal in the two
% cases, where the positive real part pole is the growth rate gamma, and p2
% is a negative real part pole. What happens when we close the loop as
% L(s) = KiGI(s) + Kz(GZ(s) is that we get a second order polynomial at the
% numerator, which may have complex conjugate roots. We hence study the
% discriminant of this polynomial, which is itself a second order function
% of the gain Ki (in the calculations below, we take Kz out of the
% parenthesis and work with L(s) = Kz(GZ(s) + Ki/Kz GI(s)).

% Impose crossing frequency and gain margin
wc = 40;
%   mg = 2;

p = roots([Gi^2 -2*(z1*Gi*Gz+2*gamma1*Gi*Gz) z1^2*Gz^2]); % Discriminant of the numerator

Fvert0 = [-p(2)*1 1]*sysZ; % we need to impose that Gi is between 0 and the root of the polynomial... (check conditions with alfredo)

% dc gain of Fvert0
% [m,~] = bode(Fvert0,0); % 
m = abs(dcgain(Fvert0));
Kz = 2.1/m; % Impone > 6db di margine di guadagno
Ki = -Kz*p(2);
FvertF = [Ki Kz]*sysZ; % Ki = -Kz*p(2)

% By removing the small damping factor zeros in the loop function
% numerator, we get rid of the knot in the Nichols diagram. Now, we are
% left with a nice system with a very high crossing frequency, which may
% negatively interact with the actual power supplies and diagnostics.
% Hence, we reduce the crossing frequency by attenuating the system at the
% desired wc (we have a very large phase margin, so the lag is not a
% problem).
%
% We use a lag network in the form (1+j*w*tau/m)/(1+j*w*tau). The network
% is parameterized in tau and m, tau being related to the characteristic
% frequency of the pole and m to the pole-zero distance. From the Nichols
% charts, we fix m = 3 (alpha below), and scan the values of tau so that we
% get the desired attenuation at the crossing frequency.

[m,~] = bode(FvertF,wc);

alpha = 3;
tau = logspace(log10(0.6/wc),log10(100/wc),50);
M   = tau*0;
phi = tau*0;
for i = 1:length(tau)
  [M(i), phi(i)] = bode(tf([tau(i)/alpha 1],[tau(i) 1]),wc);
end
[~,idx] = min(abs(M-1/m)); % find minimum distance between lag network module at wc and 1/m

R = tf([tau(idx)/alpha 1],[tau(idx) 1]); % correcting network tf
K = ss([Ki Kz]); % Full regulator

K.inputname  = {'I_{VS3}','zdot'};
K.outputname = {'V_{VS3}'};
R.inputname  = {'V_{VS3}'};
R.outputname = {'V_{VS3}'};

% sysOL = series(sysZ,K,'name');
  