 function [b,lam,Rii,yclean]=vbfa_seki(y,nl,nem);
 
% vb factor analysis
%  
% Rii: interference covariance matrix
% nk: data dimension (number of sensors)
% nt: number of time points
% nl: number of factors
% y(nk,nt): data
% b(nk:nl): mixing matrix
% lam(nk,1): sensor noise precision (diagonal elements)

% number of EM iterations
%
% initialization
%
% disp('entering reg_vbfa_seki');
 
 nk=size(y,1);  nt=size(y,2);

 Ryy=y*y';
 [U S V]=svd(Ryy/nt);
 sd=diag(S);
 b=U*diag(sqrt(sd));
 b=b(:,1:nl);
 lam=1 ./diag(Ryy/nt);
 Lam=diag(lam);
 alfa=ones(nl,1);
 iPsi=eye(nl,nl);
%
% EM iteration
%
% matrix preallocation
%
  Gamma=zeros(nl,nl);
  like=zeros(1,nem);
  
 for iem=1:nem
%
%  e-step
%
  Alfa=diag(alfa);

%   Gamma=b'*Lam*b+eye(size(Gamma));
   Gamma=b'*Lam*b+eye(size(Gamma))+nk*iPsi;
   xbar=inv(Gamma)*b'*Lam*y;
%
% calculate sufficient statistics
%
  Rxx=xbar*xbar'+nt*inv(Gamma);
  Rxy=xbar*y';
  Ryx=Rxy';
%
% m-step updates
%
  Psi=Rxx+Alfa;
  iPsi=inv(Psi);
  b=Ryx*iPsi;
  
 ilam=diag(Ryy-b*Rxy)/(nt+nl);
 lam= 1 ./ilam;
  Lam=diag(lam);
%  Lam=trace(lam)*eye(size(lam));
  
 ialfa=diag(b'*Lam*b/nk)+diag(iPsi);
 alfa=1 ./ialfa;
%
% calculating likelihood
%
% like(iem)=nt*log(det(diag(lam))/det(Gamma))-trace(diag(lam)*Ryy)+trace(Gamma*Rxx);
%
 end

% figure(1), plot(like)

% lam=trace(lam);

Rii=b*b'+diag(1./lam);

yclean=b*xbar;

% c_y=b*Rxx*b'+diag(ilam)*trace(Rxx*inv(Rxx+Alfa));
% c_y=Ryy-b*Ryx'-Ryx*b'+( b*Rxx*b'+diag(ilam)*trace(Rxx*iPsi) );
