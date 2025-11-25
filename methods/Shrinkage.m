function [sigma,alpha] = Shrinkage(R,Cor,shrink)

%See reference 4?


% Shrinkage applies the shrinkage to Cor with
% the choice of shrinkage parameter alpha.
% Sigma = alpha*F+(1-alpha)*Cor
%
% R      : The return matrix
% Cor    : The related correlation matrix
% Shrink : If it is specified alpha = shrink.

% de-mean returns
[t,n]=size(R);
meanx=mean(R);
R=R-meanx(ones(t,1),:);

% compute prior
var=diag(Cor);
sqrtvar=sqrt(var);
rBar=(sum(sum(Cor./(sqrtvar(:,ones(n,1)).*sqrtvar(:,ones(n,1))')))-n)/(n*(n-1));
F=rBar*sqrtvar(:,ones(n,1)).*sqrtvar(:,ones(n,1))';
 F(logical(eye(n)))=var;

if (nargin < 3 | shrink == -1) % compute shrinkage parameters and constant
			       
  % what we call pi-hat
  y=R.^2;
  phiMat=y'*y/t - 2*(R'*R).*Cor/t + Cor.^2;
  phi=sum(sum(phiMat));
  
  % what we call rho-hat
  term1=((R.^3)'*R)/t;
  help = R'*R/t;
  helpDiag=diag(help);
  term2=helpDiag(:,ones(n,1)).*Cor;
  term3=help.*var(:,ones(n,1));
  term4=var(:,ones(n,1)).*Cor;
  thetaMat=term1-term2-term3+term4;
  thetaMat(logical(eye(n)))=zeros(n,1);
  rho=sum(diag(phiMat))+rBar*sum(sum(((1./sqrtvar)*sqrtvar').*thetaMat));
  
  % what we call gamma-hat
  gamma = norm(Cor-F,'fro')^2;
  
  % compute shrinkage constant
  kappa=(phi-rho)/gamma;
  alpha=max(0,min(1,kappa/t));
  
else % use specified constant
  alpha=shrink;
end

% compute the estimator
sigma=alpha*F+(1-alpha)*Cor;
