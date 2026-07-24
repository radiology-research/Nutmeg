function champ_out = nut_Awsmchamp(Lp,data,flags)

% function to run Champagne using the awsm_champ 
% convention:  Y_post = Fs + Bu + v
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%inputs:
%Lp   - leadfield matrix: (sensors x dir x voxels)
%data - uses data.y for sensor data: (time x sensors x trials), data.latency
%flags - control/active time markers plus champ.* options:
%   flags.tb.timeptc  = [start stop] sample indices of the control/baseline window
%   flags.tb.timepta  = [start stop] sample indices of the active window
%   flags.cn          = 1 to column-normalize the leadfield, 0 to skip
%   flags.champ.ax     = noise-covariance mode, 0-5 (same semantics as run_champ_code.m)
%   flags.champ.multf  = scalar used by ax modes 1,2,3,4,5 (meaning depends on ax)
%   flags.champ.nem_ch = number of Champagne EM iterations
%   flags.champ.vcs    = voxel covariance structure: 0 = scalar, 1 = diagonal, 2 = general (default 2)
%   flags.champ.init   = hyperparameter init: 0 = identity+jitter (default), 1 = beamformer (inv_filter_smv1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%output
%champ_out.sources     - source time-courses (dir x voxels x time x trials)
%champ_out.sourcesmean - trial-averaged source time-courses (dir x voxels x time)
%champ_out.hyper       - voxel covariance matrices (dir x dir x voxels)
%champ_out.hyper1      - per-voxel hyperparameter (sum of diag(hyper) across dir)
%champ_out.pow         - power per voxel, averaged over the active window and trials
%champ_out.timepts     - latency vector for the active window
%champ_out.W1          - raw weight matrix from awsm_champ (interleaved dir x voxel rows)
%champ_out.wts         - de-interleaved weights (dir x sensors x voxels)
%champ_out.mn          - scale factor applied to data.y
%
%Also saves the same s_beam_CHAMP*/weights_CHAMP*/sources_CHAMP.mat files
%that run_champ_code.m produces, if a global nuts struct with
%voxels/coreg/voxelsize is available. Kept this inorder to confirm data is
%valid in similar way to clinical
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global nuts

col_norm = flags.cn;
AX       = flags.champ.ax;
multf    = flags.champ.multf;
nem_ch   = flags.champ.nem_ch;

if isfield(flags.champ,'vcs') && ~isempty(flags.champ.vcs)
    vcs = flags.champ.vcs;
else
    vcs = 2; % default: general voxel covariance 
end

if isfield(flags.champ,'init') && ~isempty(flags.champ.init)
    init = flags.champ.init;
else
    init = 0; % default: identity + jitter
end

if(AX==1||AX==2)
    display(['using ' num2str(multf) ' as scalar in noise estimate'])
elseif(AX==3)
    display(['using ' num2str(multf) ' as scaling factor for vbfa_seki0 noise estimate'])
elseif(AX==4||AX==5)
    display(['using top ' num2str(multf) ' singular vectors to build noise-estimate residual'])
end

lf = size(Lp,2); % number of orientations (2 or 3)
nv = size(Lp,3);
ns = size(data.y,2);
ntrials = size(data.y,3);

if nem_ch<50
    warning('That is not enough iterations.  We recommend more than 50.')
end

%column normalize lead field
if col_norm
    [Lp, foo] = norm_leadf_col(Lp);
end

%make data on scale 1-10
m=max(max(max(abs(data.y))));
data.y=data.y*(1/m);
data.y=double(data.y);

%extract pre- and post-stim data, per trial
nt_pre  = flags.tb.timeptc(2)-flags.tb.timeptc(1)+1;
nt_post = flags.tb.timepta(2)-flags.tb.timepta(1)+1;

prestim  = zeros(ns,nt_pre,ntrials);
poststim = zeros(ns,nt_post,ntrials);
for tr=1:ntrials
    pre  = data.y(flags.tb.timeptc(1):flags.tb.timeptc(2),:,tr)';
    pre  = pre-mean(pre,2)*ones(1,size(pre,2));
    post = data.y(flags.tb.timepta(1):flags.tb.timepta(2),:,tr)';
    post = post-mean(post,2)*ones(1,size(post,2));
    prestim(:,:,tr)  = pre;
    poststim(:,:,tr) = post;
end

timepts=data.latency(flags.tb.timepta(1):flags.tb.timepta(2));

%concatenate trials for noise-covariance estimation and Champagne
prestim2d  = reshape(prestim, ns,nt_pre*ntrials);
poststim2d = reshape(poststim,ns,nt_post*ntrials);

%per-direction leadfield slices
Lf1 = squeeze(Lp(:,1,:));
Lg1 = squeeze(Lp(:,2,:));
if lf==3
    Lh1 = squeeze(Lp(:,3,:));
else
    Lh1 = [];
end

%initialize hyperparameter matrix
if(init==0)
    ialp=zeros(lf,lf,nv);
    for i=1:nv
        for j=1:lf
            for k=1:lf
                if j==k
                    ialp(j,k,i)=1;
                else
                    ialp(j,k,i)=0;
                end
            end
        end
    end
    ialp0=zeros(lf,lf,nv);
    for iv=1:nv
        a0=randn(lf,lf);ialp0(:,:,iv)=a0*a0';
    end
    ialp1=max(max(max(ialp)));
    ialp_init=ialp+ialp0*ialp1/10000; %moves values slightly off
elseif(init==1)
    % beamformer-based init. inv_filter_smv1 only returns a 2x2 covariance
    % block per voxel. only useful for lf==2
    RzzT=poststim2d*poststim2d'/size(poststim2d,2);
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
else
    error('flags.champ.init must be 0 or 1')
end

display('computing noise covariance from pre-stim data using SEFA')
nl=20; %num of factors
nem_sefa=100;
[a1,b,lam,alp,bet,xbar1]=sefa0(poststim2d,prestim2d,nl,25,nem_sefa,0);

%compute noise covariance
if(AX==0)
    Sigma_e = double(b*b' + inv(lam));
elseif(AX==1)
    Sigma_e = multf*eye(ns,ns);
    poststim2d = a1*xbar1;
elseif(AX==2)
    Sigma_e = multf*eye(ns,ns);
elseif(AX==3)
    [b,lam2,sig,yclean] = vbfa_seki0(poststim2d,nl,nem_sefa);
    Sigma_e = multf * double(inv(diag(lam2)));
elseif(AX==4||AX==5)
    [p,d,q] = svd(poststim2d*poststim2d');
    p1 = p(:,1:multf);
    y1 = (p1*p1')*poststim2d;
    u  = poststim2d-y1;
    if(AX==4)
        Sigma_e = double(u*u');
    elseif(AX==5)
        Sigma_e = double(diag(diag(u*u')));
    end
else
    error('you somehow didn''t choose a valid AX')
end

%Leadfield
LF = zeros(ns, nv*lf);
LF(:,1:lf:end) = Lf1;
LF(:,2:lf:end) = Lg1;
if lf==3
    LF(:,3:lf:end) = Lh1;
end

display('running Champagne')
[gam,s_bar,w] = awsm_champ(poststim2d,LF,Sigma_e,nem_ch,lf,vcs,0);

%compute trial-resolved power
pow_t = 1:nt_post;
sources     = zeros(lf,nv,nt_post,ntrials);
sourcesmean = zeros(lf,nv,nt_post);
wts         = zeros(lf,ns,nv);
pow         = zeros(nv,1);

for dd=1:lf
    sd = real(s_bar(dd:lf:end,:));
    sd = sd(1:nv,:);
    sd_trials = reshape(sd,[nv nt_post ntrials]);
    sources(dd,:,:,:)    = sd_trials;
    sourcesmean(dd,:,:)  = mean(sd_trials,3);
    pow = pow + sum(sum(sd_trials(:,pow_t,:).^2,3),2)./(length(pow_t)*ntrials);
    wts(dd,:,:) = w(dd:lf:(lf*nv)-(lf-dd),:)';
end
pow = pow/lf;

%output the source time courses and hyperparameters
champ_out.sources     = sources;
champ_out.sourcesmean = sourcesmean;
champ_out.hyper       = gam;
champ_out.timepts     = timepts;
champ_out.pow         = pow;
champ_out.W1          = w;
champ_out.wts         = wts;
champ_out.mn          = m;

%get one hyperparameter value per voxel by summing across directions
temp=gam(1,1,:);
for i=2:lf
    temp=temp+gam(i,i,:);
end
champ_out.hyper1=squeeze(temp);
clear temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save
if ~isempty(nuts) && isfield(nuts,'voxels')
    voxels=nuts.voxels;
    voxelsize=nuts.voxelsize;
    coreg=nuts.coreg;
    bands=[1 160];
    if isfield(nuts,'meg') && isfield(nuts.meg,'srate')
        srate=nuts.meg.srate;
    else
        srate=1200;
    end

    timepts=1;
    sa{1}=champ_out.hyper1;
    save(['s_beam_CHAMP' num2str(vcs) 'm' num2str(multf) '_hyper.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands')

    sa{1}=real(pow);
    save(['s_beam_CHAMP' num2str(vcs) 'm' num2str(multf) '_power.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands')

    timepts=champ_out.timepts;
    sa{1}=[];
    for dd=1:lf
        sa{1}(:,:,1,dd)=squeeze(sourcesmean(dd,:,:));
    end
    save(['s_beam_CHAMP' num2str(vcs) 'm' num2str(multf) '_time.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands')

    w=wts;
    save(['weights_CHAMP' num2str(vcs) 'm' num2str(multf) '.mat'],'w')

    sources=reshape(sources,lf,nv,nt_post*ntrials);
    save('sources_CHAMP.mat','sources','-v7.3')
else
    display('No nuts.voxels found in global nuts: skipping save-to-disk of s_beam/weights/sources files.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L,colnorm] = norm_leadf_col(L)

L = permute(L,[1 3 2]);

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

L = permute(L,[1 3 2]);

%%
function [a,b,lam,alp,bet,xbar]=sefa0(y,y0,nl,n_inf,nem_init,ifplot);

% vb-em algorithm for inferring the sefa analysis model   y = a*x + b*u + v
% learn a,b,lam
%
% y(nk,nt) = data
% a(nk,nl) = mixing matrix
% lam(nk,nk) = diagonal noise precision matrix
% alp(nl,nl) = diagonal hyperparmaeter matrix
% xbar(nl,nt) = posterior means of the factors
%
% nk = number of data points
% nt = number of time points
% nl = number of factors
% nem_init = number of em iterations

nk=size(y,1);
nt=size(y,2);
nt0=size(y0,2);



disp('VBFA initialization of SEFA0');
b_init=0;
lam_init=0;
[b,lam,bet,ubar]=vbfa(y0,n_inf,nem_init,b_init,lam_init,ifplot);



a_init=0;
ryy=y*y';
ryy0=y0*y0';
if a_init==0
%   [p d q]=svd(ryy/nt);d=diag(d);
%   a=p*diag(sqrt(d));
%   a=a(:,1:nl);
   sig0=b*b'+diag(1./diag(lam));
   [p0 d0 q0]=svd(sig0);d0=diag(d0);
   s=p0*diag(sqrt(d0))*p0';
   invs=p0*diag(1./sqrt(d0))*p0';
   [p d q]=svd(invs*ryy*invs/nt);d=diag(d);
%   a=s*p(:,1:nl)*diag(sqrt(max(d(1:nl)-1,0)));
   a=s*p(:,1:nl)*diag(sqrt(abs(d(1:nl)-1)));
else
   a=a_init;
end



% initialize by svd
if 1>2
ryy=y*y';
[p d q]=svd(ryy/nt);d=diag(d);
a=p*diag(sqrt(d));
a=a(:,1:nl);

ryy0=y0*y0';
[p d q]=svd(ryy0/nt0);d=diag(d);
b=p*diag(sqrt(d));
b=b(:,1:n_inf);
lam=diag(nt0./diag(ryy0));
end



alp=diag(1./diag(a'*lam*a/nk));
alp=min(diag(alp))*diag(ones(nl,1));
bet=diag(1./diag(b'*lam*b/nk));
bet=min(diag(bet))*diag(ones(n_inf,1));

ab=[a b];
alpbet=diag([diag(alp);diag(bet)]);
nlm=nl+n_inf;
psi=eye(nlm)/(nt+nt0);

% em iteration

like=zeros(nem_init,1);
alapsi=ab'*lam*ab+nk*psi;

for iem=1:nem_init
   gam=alapsi+eye(nlm);
   igam=inv(gam);
   xubar=igam*ab'*lam*y;

   b=ab(:,nl+1:nlm);
   psib=psi(nl+1:nlm,nl+1:nlm);
   gam0=b'*lam*b+nk*psib+eye(n_inf);
   igam0=inv(gam0);
   ubar0=igam0*b'*lam*y0;

   ldlam=sum(log(diag(lam/(2*pi))));
   ldgam=sum(log(svd(gam)));
   ldgam0=sum(log(svd(gam0)));
   ldalpbet=sum(log(diag(alpbet)));
   ldpsi=sum(log(svd(psi)));
   like0=.5*nt0*(ldlam-ldgam0)-.5*sum(sum(y0.*(lam*y0)))+.5*sum(sum(ubar0.*(gam0*ubar0)));
   like(iem)=.5*nt*(ldlam-ldgam)-.5*sum(sum(y.*(lam*y)))+.5*sum(sum(xubar.*(gam*xubar)))+.5*nk*(ldalpbet+ldpsi)+like0;
if(ifplot)
   subplot(3,3,1);plot((1:iem)',like(1:iem));title('SEFA0');
   subplot(3,3,4);plot([mean(ab.^2,1)' 1./diag(alpbet)]);
   subplot(3,3,7);plot(1./diag(lam));
   drawnow;
end
   rxuxu=xubar*xubar'+nt*igam;
   ruu0=ubar0*ubar0'+nt0*igam0;
   rxuxu(nl+1:nlm,nl+1:nlm)=rxuxu(nl+1:nlm,nl+1:nlm)+ruu0;
   psi=inv(rxuxu+alpbet);

   ryxu=y*xubar';
   ryu0=y0*ubar0';
   ryxu(:,nl+1:nlm)=ryxu(:,nl+1:nlm)+ryu0;
   ab=ryxu*psi;
   lam=diag((nt+nt0)./diag(ryy+ryy0-ab*ryxu'));
   alapsi=ab'*lam*ab+nk*psi;
   alpbet=diag(nk./diag(alapsi));
end

a=ab(:,1:nl);
b=ab(:,nl+1:nlm);
alp=alpbet(1:nl,1:nl);
bet=alpbet(nl+1:nlm,nl+1:nlm);
xbar=xubar(1:nl,:);

%%
function [a,lam,alp,xbar]=vbfa(y,nl,nem_init,a_init,lam_init,ifplot);

% vb-em algorithm for inferring the factor analysis model   y = a*x + v
%
% y(nk,nt) = data
% a(nk,nl) = mixing matrix
% lam(nk,nk) = diagonal noise precision matrix
% alp(nl,nl) = diagonal hyperparmaeter matrix
% xbar(nl,nt) = posterior means of the factors
%
% nk = number of data points
% nt = number of time points
% nl = number of factors
% nem_init = number of em iterations

nk=size(y,1);
nt=size(y,2);

% initialize by svd

ryy=y*y';
if a_init==0
   [p d q]=svd(ryy/nt);d=diag(d);
   a=p*diag(sqrt(d));
   a=a(:,1:nl);
   lam=diag(nt./diag(ryy));
else
   a=a_init;
   lam=lam_init;
end

alp=diag(1./diag(a'*lam*a/nk));
alp=min(diag(alp))*diag(ones(nl,1));
psi=eye(nl)/nt;

% em iteration

like=zeros(nem_init,1);
alapsi=a'*lam*a+nk*psi;


for iem=1:nem_init
   gam=alapsi+eye(nl);
   igam=inv(gam);
   xbar=igam*a'*lam*y;


if(ifplot)
   subplot(3,3,1);plot((1:iem)',like(1:iem)/nt);title('VBFA: like');
   subplot(3,3,4);plot([mean(a.^2,1)' 1./diag(alp)]);title('1/alp');
   subplot(3,3,7);plot(1./diag(lam));title('1/lam');
   drawnow;
end

   ryx=y*xbar';
   rxx=xbar*xbar'+nt*igam;
   psi=inv(rxx+alp);

   a=ryx*psi;
   lam=diag(nt./diag(ryy-a*ryx'));
   alapsi=a'*lam*a+nk*psi;
   alp=diag(nk./diag(alapsi));
end

return
%%
