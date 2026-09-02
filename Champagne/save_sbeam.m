function save_sbeam(gamma,timecourse,t1,t2,nuts)

%inputs

%gamma: the hyperparametrs (numvox x 1)
%timecourse: the time courses (numvox x numtime x numdirs)
%t1: start of post-stim window
%t2: stop of post-stim window
%nuts: nuts from session file

coreg=nuts.coreg;
srate = nuts.meg.srate;
voxels=nuts.voxels;
voxelsize = nuts.voxelsize;
bands = [1 100];
timepts=1;

sa{1} = gamma;

save(['s_beam_hyper.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands')

sa=[];

tt1=dsearchn(nuts.meg.latency,t1);
tt2=dsearchn(nuts.meg.latency,t2);
timepts=[];
timepts = nuts.meg.latency(tt1:tt2);


sa{1}(:,:,1,1)=timecourse(:,:,1);
if size(timecourse,3)==2
  sa{1}(:,:,1,2)=timecourse(:,:,2);  
elseif size(timecourse,3)==3
sa{1}(:,:,1,3)=timecourse(:,:,3);
end

save(['s_beam_time.mat'],'sa','coreg','srate','timepts','voxels','voxelsize','bands')

