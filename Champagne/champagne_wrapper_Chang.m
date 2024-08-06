
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Changpang_wrapper for epilepsy spike or aef or SEF localization by Chang Cai
%
% 2020.01.09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
% save the result in path_results
%
% Input:
% path_session: the path and name of the session file
% start_tp: the start time point of the spike
% end_tp: the end time of the spike
% trial_idx: which trial
% avg_flag: to average the trials or not. flag=1:avg;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function champagne_wrapper_Chang(path_session,start_tp,end_tp,trial_idx,avg_flag)

if nargin < 5, error('Please set the right parameters'); end

%%
thisdir = pwd;
load(path_session);

%% parameters setting
[nc nd nv] = size(Lp);
nem = 100;
vcs = 0;
plot_on =1;
coup = 0;
ncf = 1;
t1 = find(meg.latency==start_tp);
t2 = find(meg.latency==end_tp);
timepts = meg.latency(t1:t2);
bands = [1 160];  % default setting
srate = meg.srate;

%% leadfield normalization
f = reshape(Lp,nc,nd*nv);
for i = 1:nd*nv
    f(:,i) = f(:,i)./sqrt(sum(f(:,i).^2));
end

%% data setting
y = meg.data(:,:,trial_idx);

[np nc nt]= size(y);
for i=1:nc
    for j=1:nt
        y(:,i,j) = detrend(y(:,i,j));
    end
end

%% output path
[pathstr,sessname] = fileparts(path_session);
output_dir = fullfile(pathstr,sprintf('%s_post%dto%dms',sessname,start_tp,end_tp));
bSave = savecheck(output_dir,'dir');

if ~bSave
    disp('Directory exists: Champagne not run.')
    return
end
mkdir(output_dir);

if avg_flag == 1
    
    ypost = mean(y(t1:t2,:,:),3)';
    
    %% running Champagne with noise learning
    [gamma,x,w,c,l, sigu,v]=awsm_champ_noiseup(ypost,f,nem,nd,vcs,plot_on,coup,ncf);
    
    
    %% save results which can be viewed through nutmeg.
    sa = [];
    sa{1}(:,:,1,1) = x(1:nd:end,:);
    if nd>1
        sa{1}(:,:,1,2) = x(2:nd:end,:);
    end
    if nd>2
        sa{1}(:,:,1,3) = x(3:nd:end,:);
    end
    
    cd(output_dir);
    %     resultname = 's_beam_CHAMP_noiselearningavg.mat';
    %     %%
    %     bSave = savecheck(output_dir,'dir');
    %
    %     if ~bSave
    %         disp('Directory exists: Champagne not run.')
    %         return
    %     end
    %%
    
    save(['s_beam_CHAMP_noiselearningavg.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands')
    cd(thisdir);
else
    for i_trial=1:nt
        ypost = squeeze(y(t1:t2,:,i_trial))';
        
        %% running Champagne with noise learning
        [gamma,x,w,c,l, sigu,v]=awsm_champ_noiseup(ypost,f,nem,nd,vcs,plot_on,coup,ncf);
        
        %% save results which can be viewed through nutmeg.
        sa = [];
        sa{1}(:,:,1,1) = x(1:nd:end,:);
        if nd>1
            sa{1}(:,:,1,2) = x(2:nd:end,:);
        end
        if nd>2
            sa{1}(:,:,1,3) = x(3:nd:end,:);
        end
        
        cd(output_dir);
        save(['s_beam_CHAMP_noiselearning_trial',num2str(trial_idx(i_trial)),'.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands')
        cd(thisdir);
    end
    
end



end