

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learn both voxel covariances and noise covariance in Champagne using
% augmentd leadfield
% 
% � 2019 Convex Imaging by Chang Cai
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Output: 
% gamma(nd,nd,nv) = voxel covariance matrices
% x(nv*nd,nt) = voxel current density 
% w(nv*nd,nk) = reconstruction filter
%
% Input:
% y(nk,nt) = sensor data 
% f(nk,nv*nd) = lead field matrix
% nd = # of orientations
% vcs = voxel covariance structure: 0 = scalar, 1 = diagonal
% nem = maximum # of iterations
%
%% add two more parameter for input
% coup: covaraince update methods, 0=covex boounding, 1 = EM, 2 = Mackay
% ncf: noise covavariance form, 0 = scalar, 1 = heter
%
% #s:
% nk = # of sensors
% nv = # of voxels
% nt = # of data time points




function [gamma,x,w,c,l, sigu,v]=awsm_champ_noiseup(y,f,nem,nd,vcs,plot_on,coup,ncf);

[gamma x w l v c,sigu]=champ_vec(y,f,nem,nd,vcs,plot_on,coup,ncf);

return




function [gamma,x,w,like,vvec,c,sigu]=champ_vec(y,f,nem,nd,vcs,plot_on,coup,ncf);


eps1=1e-8;
[nk nvd]=size(f);
nv=nvd/nd;
nt=size(y,2);

%% augmented leadfield 
nv_old = nv;
nvd_old = nvd;

lf_n = reshape(repmat(eye(nk),nd,1),nk,nk*nd);
f = [f lf_n];

[nk nvd]=size(f);
nv=nvd/nd;

cyy=y*y'/nt;

% Initialize voxel variances

f2=sum(f.^2,1);             %which filter used to estimated 
invf2=zeros(1,nvd);
ff=find(f2>0);
invf2(ff)=1./f2(ff);
f=double(f);
w=spdiags(invf2',0,nvd,nvd)*f';
%w=spdiags(1./sum(f.^2,1)',0,nvd,nvd)*f';
inu0=mean(mean((w*y).^2));
%gamma=inu0*repmat(eye(nd,nd),[1,1,nv]);
vvec=inu0*ones(nvd,1);


% Learn voxel variances
if(plot_on)
    Champfig=uifigure('AutoResizeChildren','off');
end

like=zeros(nem,1);
for iem=1:nem
%    iem
    vmat=spdiags(vvec,0,nvd,nvd);
    c=f*vmat*f';   %used for 4.74
%    [p d q]=svd(c);
    [p d]=eig(c);%used for 4.74
    d=max(real(diag(d)),0);%used for 4.74
    invd=zeros(nk,1);%used for 4.74
    ff=find(d>=eps1);%used for 4.74
    invd(ff)=1./d(ff);   %concave function parameters 1/a1,1/a2,....1/aN %used for 4.74
    invc=p*spdiags(invd,0,nk,nk)*p';        %model data convariance y %used for 4.74
    %invc = inv(c);
 
    like(iem)=-.5*(sum(log(max(d,eps1)))+nk*log(2*pi))-.5*sum(sum(invc.*cyy));  %4.77 4.18 y's distribution directly

    if(plot_on)
        p1=subplot(2,2,1,'Parent',Champfig);
        plot(p1,(1:iem),like(1:iem));
        title(['Likelihood: ' int2str(iem) ' / ' int2str(nem)]);
        xlabel('iteration');
        set(p1,'XLim',[0 iem]);
    end
    
    fc=f'*invc;         %4.82
    w=vmat*fc;          %4.82 update s
    x=w*y;              %4.82 source distribution
    x2=mean(x.^2,2);    %xk distribution is equal to the posterior mean of the source distribution
    z=sum(fc.*f',2);    %update Z 4.79
    
    if 0 == coup  %% convex update
        if vcs==0
            x20=sum(reshape(x2,nd,nv),1);
            z0=sum(reshape(z,nd,nv),1);
            v0=zeros(size(z0));
            ff=find(z0>0);
            v0(ff)=sqrt(max(x20(ff)./z0(ff),0)); % CN 10/2012 added max,0
    %        v0(ff)=sqrt(x20(ff)./z0(ff));
            vvec=reshape(ones(nd,1)*v0,nvd,1);
        else
            vvec=zeros(size(x2));
            ff=find(z>0);
            vvec(ff)=sqrt(max(x2(ff)./z(ff),0)); % CN 10/2012 power 4.88 update v
        end
        vvec = real(vvec);
        v=sum(reshape(vvec,nd,nv),1);   %voxel power
    elseif 1 == coup %% EM update
        
%         ff = find(vvec>0);
%         igam = zeros(size(vvec));
%         igam(ff) = 1./vvec(ff);
%         
        igam = sum((eye(nvd)-w*f).*vmat,2);
        vvec = (igam)+x2;
        v = sum(reshape(vvec,nd,nv),1);   %voxel power
    elseif 2 == coup %% Mackay
        ZZ = sum(f'*invc.*f',2);%insigu*invcovx*f'*f;
        dzz = vvec.*ZZ;
        ff = find(dzz>0);
        vvec(ff) = x2(ff)./dzz(ff);
        v = sum(reshape(vvec,nd,nv),1);   %voxel power
    end
    
        
    if 0==vcs
        vvec = mean(reshape(vvec,nd,nv),1);
        vvec = reshape(repmat(vvec,nd,1),nvd,1);
    end
    
    if 0==ncf
        vvec(nvd_old+1:end) = mean(vvec(nvd_old+1:end));
        v = sum(reshape(vvec,nd,nv),1);   
    end
    sigu = diag(vvec(nvd_old+1:end));
    
    if(plot_on)
        p2=subplot(2,2,2,'Parent',Champfig);
        plot(p2,(1:nv),v);
        title(['Voxel power: ' num2str(nv) ' / ' num2str(nv)]);
        xlabel('voxel index');
        set(p2,'XLim',[1 nv]);
        drawnow
    end
    
end
   

sigu = diag(mean(reshape(vvec(nvd_old+1:end),nd,nk),1));
vvec = vvec(1:nvd_old);
x = x(1:nvd_old,:);
nv = nv_old;
w = w(1:nvd_old,:);

if nd==1
    gamma=reshape(vvec,1,1,nv);
else
    gamma=zeros(nd,nd,nv);
    for iv=1:nv
        gamma(:,:,iv)=diag(vvec((iv-1)*nd+1:iv*nd));
    end
end

return


