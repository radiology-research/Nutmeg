function [s_hat,Smat]=inv_filter_smv1(InvRzz,Lf1,Lg1,Lh1)
%
%  implement scalar minimum-variance spatial filter
%  for initializing Champagne
%  
%  Kensuke Sekihara  July 31, 2010
%
%  Nc: number of channels
%  Np: number of pixels
%  Ne: number of trials (epochs)
%
% inputs:
% Lf1, Lg1: lead field    LKf1(1:Nc,1:Np)
% InvRzz: inverse covariance matrix    InvRzz(1:Nc,1:Nc)
%
% outputs:
% s_hat: spatial -filter power outputs , S(1:Np,2)
% Smat: source power matrix Smat(2,2,1:Np)
% 
  Np=size(Lf1,2);
%
%  matrix preallocation ------------------------
%
 s_hat=zeros(Np,1);
 Smat=zeros(2,2,Np);
%
%h = waitbar(0,'spatial-filter scanning');
%          
  for nvx=1:Np
%
        lpf=Lf1(:,nvx); 
        lpg=Lg1(:,nvx);
        if(~isempty(Lh1))
            lph=Lh1(:,nvx);
            Lp=[lpf,lpg,lph]; 
        else
        Lp=[lpf,lpg];
        end
%
%    orientation estimation
%
	 L1=Lp'*InvRzz*Lp;    
     L0=Lp'*Lp; 
     [Up,Sp]=eig(L1,L0);
     [dia_fix,iorder]=sort(abs(diag(Sp)));                
%      Up(:,:)=Up(:,iorder);
%      normalpm=Up(:,1);
%      lvp=Lp*normalpm;  lvp=lvp/norm(lvp);
% 
% source reconstruction
% 
    spow=1/dia_fix(1);
     s_hat(nvx)=spow;
%
    Smat(1,1,nvx)=spow/2;
    Smat(2,2,nvx)=spow/2;
    
%    waitbar(nvx/Np,h)
end
%       close(h)
%
