clc
clear
close all
path_sess = append(pwd,'\session_Llip.mat');

start_tp = -25;
end_tp = 125;
trial_indx = [1]; %include range of trials or 1 for averaged dataset
avg_flag = 1; % 1 means avg

champagne_wrapper_Chang(path_sess,start_tp,end_tp,trial_indx,avg_flag);
