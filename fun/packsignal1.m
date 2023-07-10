function U = packsignal1(k,t,N,u)
% Pack signal in a 2D Hankel matrix as in the paper by Persis&Tesi
[~, L] = size(u);
if k+1+t+N-2<=L
  U = zeros(numel(t),N);
  for i=(k+1):t+k
    U(i-k,:) = u(:,i:i+N-1);
  end
else
  U = [];
end
