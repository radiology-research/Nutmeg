
function [ypost0,Nee_op] = MYptsss(ypost0,Fc,mu,Nee)
% ypost0 = MYptsss(ypost0,Fc,mu,Nee)
% poorman's tSSS
%
% Nc: number of sensors
% Nt: number of time points
% Np: number of pixels
%
% ypost0(1:Nc,1:Nt): input (spatio-temporal) data
% Fc(1:Nc,1:2*Np): lead field matrix. 
% Fc(1:Nc,1:3*Np): If you use a three component lead field
%
% mu: Dimension of Lead Field Span
% Nee: Dimension of Intersection

%%
 G=Fc*Fc';           % computing gram matrix
 [U,S]=eig(G);
eig_els=abs(diag(S));
[eig_el, iorder]=sort(-eig_els);%small to large
eig_el=-eig_el;
U(:,:)=U(:,iorder);
ev_max=length(eig_el);
% figure()
% plot(1:ev_max,eig_el,'-o')
% axis([1,ev_max, 0, eig_el(1)*1.1])
% title('LF eigen value spectrum')
% xlabel('order of eigenvalues')
% ylabel('relative magnitude')


%
% mu=input('enter dimension of lead field span ? (30(E)) : ');
% close;
% if isempty(mu)==1, mu=30; end
Us = U(:,1:mu);             % signal subspace define  
USUS=Us*Us';
%
% projection onto the inside of the span
%
yin=USUS*ypost0;
%
% projection onto the outside of the span
%
yout=(eye(size(USUS))-USUS)*ypost0;

Nec=40;
Net=40;
[Ae,Nee,Nee_op]=MY_CSP0(yin,yout,Nec,Net,Nee);
%  ypost0=ypost0*(eye(size(AeAe))-AeAe);
 
 ypost0 = ypost0-(ypost0*Ae)*Ae';

