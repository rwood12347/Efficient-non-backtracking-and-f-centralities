function newMat = NewnodeNBTW(M, a, n, t)
% testinv = shortcircleinv(M,n,t,a);
  E = {};
  e = ones(n,n);
 for i = 1:t
     E{i} = e;
 end
 E = blkdiag(E{:});
testinv =circleinv(E - newcircle(trans(a*M,n,t),a*M,n,t),n,t);

 
X = full(a*newcircle(M,testinv,n,t));
test = eye(n*t) - X + dd(a*M*X,n,t);
Y = test \ eye(n*t);
newMat = Y;




end