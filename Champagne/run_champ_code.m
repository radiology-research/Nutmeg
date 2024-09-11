function [sources, ialp, gamma, pow, wts] = run_champ_code(nuts,numtrials,pre_start,pre_stop,post_start,post_stop,init,vcs,AX,multf)

iters=100;

lf=size(nuts.Lp,2);
nv=size(nuts.Lp,3);
nk = size(nuts.Lp,1);

time=nuts.meg.latency;

t1=dsearchn(time,pre_start)
t2=dsearchn(time,pre_stop)
t3=dsearchn(time,post_start)
t4=dsearchn(time,post_stop)

pre=nuts.meg.data(t1:t2,:,numtrials);
post=nuts.meg.data(t3:t4,:,numtrials);

pre=permute(pre,[2 1 3]);
for tr=1:size(pre,3)
    pre(:,:,tr)=remove_dcoffset(pre(:,:,tr));
end
%cMEGdata=1e14*detrend(reshape(pre,size(pre,1),size(pre,2)*size(pre,3)));
cMEGdata=detrend(reshape(pre,size(pre,1),size(pre,2)*size(pre,3))); %CN 3/2012


post=permute(post,[2 1 3]);
for tr=1:size(post,3)
    post(:,:,tr)=remove_dcoffset(post(:,:,tr));
end
%MEGdata=1e14*detrend(reshape(post,size(post,1),size(post,2)*size(post,3)));
MEGdata=detrend(reshape(post,size(post,1),size(post,2)*size(post,3))); %CN 3/2012

if(lf==2)
    Lf1=squeeze(nuts.Lp(:,1,:));
    Lg1=squeeze(nuts.Lp(:,2,:));
    [Lf1,colnorm] = norm_leadf_col(Lf1);
    [Lg1,colnorm] = norm_leadf_col(Lg1);
    Lh1=[];
else
    Lf1=squeeze(nuts.Lp(:,1,:));
    Lg1=squeeze(nuts.Lp(:,2,:));
    Lh1=squeeze(nuts.Lp(:,3,:)); %CN 8/2012 (was: (:,2,:)

    [Lf1,colnorm] = norm_leadf_col(Lf1);
    [Lg1,colnorm] = norm_leadf_col(Lg1);
    [Lh1,colnorm] = norm_leadf_col(Lh1);
end


if(init==0)
    for i=1:nv
        for k=1:lf
            for j=1:lf
                if k==j
                    ialp(k,j,i)=1;
                else
                    ialp(k,j,i)=0;
                end
            end
        end
    end
    %puts random # in iapl0 with off diags the same
    ialp0=zeros(lf,lf,nv);
    for iv=1:nv
        a0=randn(lf,lf);ialp0(:,:,iv)=a0*a0';
    end
    ialp1=max(max(max(ialp)));
    ialp_init=ialp+ialp0*ialp1/10000; %moves values sligtly off

elseif(init==1)
    RzzT=MEGdata*MEGdata'/size(MEGdata,2);
    condRzz=cond(RzzT);
    fprintf('condition number of RzzT is %.4g \n',condRzz);
    Nc=size(Lf1,1);
    rgamma=1e-4;
    if condRzz>1e6
        RzzT=RzzT+rgamma*max(eig(RzzT))*eye(Nc,Nc);
    end
    InvRzz=inv(RzzT);
    [s_hat,Smat]=inv_filter_smv1(InvRzz,Lf1,Lg1,Lh1);
    maxabs_ialp=max(max(max(abs(Smat))));
    ialp_init=Smat/maxabs_ialp;
end

ypre=cMEGdata; %CN 8/2014: ypre=double(cMEGdata);
ypost=MEGdata; %CN 8/2014: ypost=double(MEGdata);
clear cMEGdata MEGdata

tic;
nl=20;
nem=100;
[a1,b,lam,alp,bet,xbar1]=sefa0(ypost,ypre,nl,25,nem,0);

H=b;
ni=size(H,2);
lambda_init = ones(ni,1);
Sigma_f=diag(1./lam);


if(AX==0)
    Sigma_e = double(b*b' + inv(lam));
elseif(AX==1)
    %multf=double(multf*min(eig(ypost'*ypost)));
    Sigma_e = multf*eye(nk,nk);
    ypost=a1*xbar1;
elseif(AX==2)
    Sigma_e = multf * eye(nk,nk);
elseif(AX==3)
    [b,lam2,sig,yclean]=vbfa_seki0(ypost,nl,nem);
    Sigma_e = multf * double(inv(diag(lam2)));
elseif(AX==4|AX==5)
    [p d q]=svd(ypost*ypost');
    p1=p(  :,1:multf);
    y1=(p1*p1')*ypost;
    u=ypost-y1;
    if(AX==4)
        Sigma_e=double(u*u');
    elseif(AX==5)
        Sigma_e=double(diag(diag(u*u')));
    end
end

LF=zeros(size(Lf1,1),size(Lf1,2)*lf);
LF(:,1:lf:end)=Lf1;
LF(:,2:lf:end)=Lg1;
if(lf==3)
  LF(:,3:lf:end)=Lh1;  
end

tic
%LF = double(LF); %CN 08/2014 change precision

% changed plot_on from on to off (1 to 0)
[gamma,s_bar,w]=awsm_champ(ypost,LF,Sigma_e,iters,lf,vcs,0);

numtrials=length(numtrials);
t=time(t3:t4);
pt1=dsearchn(t,t3);
pt2=dsearchn(t,t4);
pow_t=[pt1:pt2];
tr=numtrials;

nt=size(s_bar,2);

ialp=squeeze(gamma(1,1,1:nv)+gamma(2,2,1:nv));

s1=s_bar(1:lf:end,:);
s11=s1(1:nv,:);
s111=reshape(s11,[nv nt/numtrials numtrials]);
pw1=sum(sum(s111(:,pow_t,1:tr).^2,3),2)./(length(pow_t)*tr);
s1mean=mean(s111,3);

s2=s_bar(2:lf:end,:);
s22=s2(1:nv,:);
s222=reshape(s22,[nv nt/numtrials numtrials]);
pw2=sum(sum(s222(:,pow_t,1:tr).^2,3),2)./(length(pow_t)*tr);
s2mean=mean(s222,3);

if(lf==3)
 s3=s_bar(3:lf:end,:);
s33=s2(1:nv,:);
s333=reshape(s33,[nv nt/numtrials numtrials]);
pw3=sum(sum(s333(:,pow_t,1:tr).^2,3),2)./(length(pow_t)*tr);
s3mean=mean(s333,3);
pow=.5*(pw1+pw2+pw3);  

sources(1,:,:)=s1;
sources(2,:,:)=s2;
sources(3,:,:)=s3;
else
    sources(1,:,:)=s1;
sources(2,:,:)=s2;
  pow=.5*(pw1+pw2);  
end

wts(1,:,:)=w(1:lf:(lf*nv)-(lf-1),:)';
wts(2,:,:)=w(2:lf:(lf*nv)-(lf-2),:)';
if(lf==3)
    wts(3,:,:)=w(3:lf:(lf*nv),:)';
end

if(isfield(nuts,'voxels'))
timepts=1;
voxels=nuts.voxels;
voxelsize=nuts.voxelsize;
coreg=nuts.coreg;
bands=[1 160];
srate=1200;
sa{1}=ialp;

save(['s_beam_CHAMP' num2str(vcs) 'm' num2str(multf) '_hyper.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands')

sa{1}=pow;

save(['s_beam_CHAMP' num2str(vcs) 'm' num2str(multf) '_power.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands')

timepts=t;
sa{1}=[];
sa{1}(:,:,1,1)=s1mean;
sa{1}(:,:,1,2)=s2mean;
if(lf==3)
    sa{1}(:,:,1,3)=s3mean;
end
save(['s_beam_CHAMP' num2str(vcs) 'm' num2str(multf) '_time.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands')

w=wts;
save(['weights_CHAMP' num2str(vcs) 'm' num2str(multf) '.mat'],'w')

save('sources_CHAMP.mat','sources','-v7.3') %CN 3/2012: added -v7.3 flag

end

function X=remove_dcoffset(X)
%input X [number sensors x number of voxels]

m=mean(X,2);

mm=m*ones(1,size(X,2));

X=X-mm;

return

function [L,colnorm] = norm_leadf_col(L)

if size(L,3)>1
    for i=1:size(L,3)
        colnorm(:,i) = sqrt(sum(L(:,:,i).^2));
        L(:,:,i) = L(:,:,i)./repmat(colnorm(:,i)',[size(L,1) 1]);
    end
else
    for i=1:size(L,2)
        colnorm(:,i) = sqrt(sum(L(:,i).^2));
        L(:,i) = L(:,i)./repmat(colnorm(:,i)',[size(L,1) 1]);
    end
end

