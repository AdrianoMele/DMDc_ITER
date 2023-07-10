function U = packsignal2(k,T,u)
% Pack signal in a 1D Hankel matrix as in the paper by Persis&Tesi

% [m, L]=size(u);
%U=zeros(m*t,N);
U=u(:,k+1:T+k);
    

