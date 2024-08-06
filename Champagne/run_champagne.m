function [sources, ialp, gamma, pow, wts] = run_champagne(directory,subj,sess_name,trials,t1,t2,t3,t4,init,vcs,ax)

%%function to run new awsm_champ.m Champagne code

%%INPUTS can either be 1) a session file from NUTMEG or 2) a structure with the
%%data, lead field, and time vector (need to change sess_name input to
%%reflect which type of input you requier)

% directory - the full path to the directory where you subject's directory 
    %is located

% subj - the name of the directory with your session file or where you want 
    %to save the output

% sess_name - either 1) the name of your NUTMEG session file or 2) a structure
% with data, lead field, and time vector

    %if using option 2) then you need to input sess_name as a struture with:
    % sess_name.data - sensor data, [num time x num sensors x num trials]
    % sess_name.LF = lead field, [num sensors x num orientations x num voxels]
    
% trials - number of trials to use (can be 1 if averaged data)

% t1 - start of pre-stim pd
% t2 - end of pre-stim pd
% t3 - start of post-stim pd
% t4 - end of post-stim pd

% init - initialization of hyperparameters, 1 = initialize with ones, 2 =
    % initialize with beamforming
    
% vcs - voxel covariance structure, 2 = full voxel covariance, 1 = diagonal 
    %voxel covariance, 0 = scalar voxel covariance
    
% ax - option to change the noise covariance (advanced setting), default
% should be [0 0]
    
%%OUTPUTS

%sources - source time courses
%ialp - sum of hyperparamteres (if applicable)
%gamma - full hyperparameter matrix
%pow - power at each voxel
%wts - weight matrix used to create sources

%example usage: [sources, ialp, gamma, pow, wts] = run_champagne('/data/research_meg8/Champagne_for_qb3/champagne/compiled_code','new_champagnecode','JPO_AEF_sess_8mm_st.mat',10,-95,-5,5,200,0,2,[0 0]);


cd(directory)
cd(subj)

if(isstruct(sess_name))
    data_all=sess_name;
  nuts.meg.data=data_all.data;  
  nuts.Lp=data_all.LF;
  nuts.meg.latency=data_all.latency;
else   
load(sess_name);
mkdir([sess_name(1:end-4) '_' num2str(length(trials)) 'trials_' num2str(t3) 'to' num2str(t4) '_AX' num2str(ax(1))])
cd([sess_name(1:end-4) '_' num2str(length(trials)) 'trials_' num2str(t3) 'to' num2str(t4) '_AX' num2str(ax(1))])
end

if(length(trials)==1)
   trials=[1:trials];
end

[sources, ialp, gamma, pow, wts] = run_champ_code(nuts,trials,t1,t2,t3,t4,init,vcs,ax(1),ax(2));
