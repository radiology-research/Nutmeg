function convert_to_old_sbeam(sbeam)

load([sbeam '.mat'])

beam.timewindow=timepts;
beam.voxels=voxels;
beam.voxelsize=voxelsize;
beam.coreg=coreg;
beam.srate=srate;
beam.s_z=0;
beam.s_th=sa{1};
beam.s_ph=0;
beam.phi=0;
beam.W1_th=0;
beam.W1_ph=0;

save([sbeam 'old.mat'],'beam')