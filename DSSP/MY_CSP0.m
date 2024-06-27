
function [Ae,Nee,Nee_op]=MY_CSP0(ysig,yint,Nec,Net,Nee)
%
% interference rejection by removing the common temporal subspace of the two subspaces
% K. Sekihara
% March 28, 2012
% Golub and Van Loan, Matrix computations, The Johns Hopkins University Press, 1996
%
%  Nc: number of channels
%  Nt: number of time points
%  yint(1:Nc,1:Nt): interference data
%  ysig(1:Nc,1:Nt): signal plus interference data
%  ypost(1:Nc,1:Nt): denoised data
%  Nec: dimension of the interference subspace
%  Net: dimension of the signal plus interference subspace
%  Nee: dimension of the intersection of the two subspaces

[~,Sc,Vc]=svd(yint,'econ');
% dd=diag(Sc);
% figure(),plot(1:50,dd(1:50),'-o');pause; close
% Nec=input('enter IFS dimension: ');
%
[~,St,Vt]=svd(ysig,'econ');
% dd=diag(St);
% figure(),plot(1:50,dd(1:50),'-o');pause; close
% Net=input('enter SIFS dimension: ');
% 
Qc=Vc(:,1:Nec);
Qt=Vt(:,1:Net);

C=Qt'*Qc;
[U,S,V]=svd(C);

%figure(); plot(diag(S),'-o')
%axis([0 50 0.9 1])
% Nee=input('enter dimension of the intersection (recommend to underestimate): ');
% close;
Nee_op = find(diag(S)>0.95,1,'last') % check for optimal Nee based on alogorithm
A1=Qt*U;
% A2=Qc*V;
if Nee_op<Nee
    Nee = Nee_op;
end
if Nee<5 | isempty(Nee_op)
    Nee =5;
end

%Nee = Nee_op;
Ae=A1(:,1:Nee);
% AeAe=Ae*Ae';
%ypost=ysig*(eye(size(AeAe))-AeAe);
Nee_op = Nee;
return