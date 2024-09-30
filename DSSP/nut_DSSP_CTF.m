function nut_DSSP_ctf(session_name, dsfoldername)
% for CTF MEG data only, will write out a cleaned session file or CTF ds
% session_name is the sessionfile name with leadfield ** required
% dsfoldername is the .ds name
%
% if only a session file (string, path/to/filename) is entered, only a cleaned sessionfile will output 
% if a sessionfile and ds dir entered (strings, path/to/target), a cleaned ds dir will output
% if no input variables, it will ask you about your data
%
% Corby L. Dale, wrapping code by Chang Cai, which originated from Kensuke Sekihara 
% this version created and tested August 1, 2024

% add UCSF RHEL9/RHEL7 paths
if isfolder('/netopt/share/bin/local/bil/')==1   % test that you are on the UCSF system
    addpath(genpath('/netopt/share/bin/local/bil/spm8_nutmeg_lhinkley/'))
    addpath(genpath('/data/research_meg/rendering_tools'))
    addpath(genpath('/netopt/share/lib/local/bil/matlab/Nutmeg-Next/'))
    addpath(genpath('/netopt/share/lib/local/bil/matlab/Nutmeg-Next/DSSP/'))
    addpath('/data/research_meg11/ccai/DSSP_Champ_ctf/ctf_fieldtrip');
    addpath('/data/research_meg11/ccai/DSSP_Champ_ctf/DSSP_new');
    addpath('/data/research_meg11/ccai/DSSP_Champ_ctf/ctf_code');
    addpath('/data/research_meg11/ccai/DSSP_Champ_ctf');
    addpath(genpath('/data/research_meg10/ccai/ccai/public_data/fieldtrip/fieldtrip-20160222'));
else
    disp('make sure your path to nutmeg is added before running this function')
end

ctftools=which('writeCTFds');
if isempty(ctftools)
    error('Path to CTF I/O tools is not set correctly')
end

narginchk(0,2); % don't bother if not called properly

thisfolder=pwd;
targfolder=thisfolder;
outputflag = 2;  % we assume a 1 argument function has a session file

if nargin <1  % then specify one or two files for input
    disp('Please specify a sessionfile.')
    [session_name, targfolder]=uigetfile('*.mat', 'Select session file');
    if session_name==0
        error('there must be a session file with a lead field to run DSSP')
    end
    disp('Please select a .ds folder or cancel to return a cleaned session file.')
    dsfoldername=uigetdir(thisfolder, 'Select .ds folder or cancel');
    if dsfoldername==0
        outputflag = 2;
    else 
        outputflag = 1;
    end
end
if nargin == 2  % default to ds cleaning if both present
    outputflag = 1;
end

mu = 50;
Nee = 10;

if outputflag == 1
    % read .ds file
    ds = readCTFds(dsfoldername);
    % get MEG channel index
    chan=cellstr(ds.res4.chanNames);
    j=0;
    for i=1:size(chan)
        if strcmp('M',chan{i}(1))
            j = j+1;
            ind(j)=i;
        end
    end
    % get CTF data
    dat_in = getCTFdata(ds); %
    dat_in_meg = dat_in(:,ind,:);
end

% load leadfield information (and MEG data if sessionfile-only)
load(fullfile(targfolder, session_name));

if outputflag == 2
    dat_in_meg = meg.data;  
end

%remove dc offset
for i=1:size(dat_in_meg,2)
    for j = 1:size(dat_in_meg,3)
        dat_in_meg(:,i,j) = detrend(dat_in_meg(:,i,j),'linear');%-mean(dat_in_meg(:,i,j));
    end
end

% prepare file to log optimal Nee for future reference
[~,name,~]=fileparts(fullfile(targfolder,session_name));
fid = fopen(strcat(thisfolder,'/', name,'_opNee.txt'),'w+');
fprintf(fid, 'This file logs the optimal Nee determined by its math definition in algorithm for furture reference.\n');

% run DSSP
data_dssp = [];
Fc = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));
[nt,nk,ntrial]=size(dat_in_meg);    

%{
fprintf(fid, 'No sample length assessment: Using whole trial for DSSP treatment\n');
for i=1:ntrial
    str = [' trial = ',num2str(i)];
    disp(str);
    [tmp,Nee_op] = MYptsss(dat_in_meg(:,:,i)',Fc, mu, Nee);
    data_dssp(:,:,i) = tmp';
    fprintf(fid,'trial %d, Nee_op = %d\n',i,Nee_op);
    close all;
end
%}

% Take sample length into account
if nt<10000
    fprintf(fid, 'Short trials, using whole trial for DSSP treatment\n');
     for i=1:ntrial
         str = [' trial = ',num2str(i)];
         disp(str);
         [tmp,Nee_op] = MYptsss(dat_in_meg(:,:,i)',Fc, mu, Nee);
         data_dssp(:,:,i) = tmp';
         fprintf(fid,'trial %d, Nee_op = %d\n',i,Nee_op);
         close all;
     end
else
     temp = round(size(dat_in_meg,1)/6000);
     fprintf(fid, 'Long trials, using 6000 data samples per DSSP treatment\n');
     for i=1:ntrial
     str = [' trial = ',num2str(i)];
     disp(str);
         for j= 1:temp-1
                 [tmp,Nee_op] = MYptsss(dat_in_meg((j-1)*6000+1:j*6000,:,i)',Fc, mu, Nee);
                  data_dssp((j-1)*6000+1:j*6000,:,i) = tmp';
                 fprintf(fid,'trial %d, Nee_op = %d\n',i,Nee_op);
         end
         [tmp,Nee_op] = MYptsss(dat_in_meg((temp-1)*6000+1:nt,:,i)',Fc, mu, Nee);
          data_dssp((temp-1)*6000+1:nt,:,i) = tmp';
         fprintf(fid,'trial %d, Nee_op = %d\n',i,Nee_op);
         close all
     end
end
%}
Nee = Nee_op;
fclose(fid); % close optimal Nee logfile

if outputflag==1  %% save CTF ds out
    dat_out = dat_in;
    dat_out(:,ind,:)=data_dssp(:,:,:);  % for CTF data (has ref chans, but do not 3grad correct denoised data!) 
    name = ds.baseName;
    ds_temp=writeCTFds(strcat(fullfile(thisfolder,name),'_Dssp','_mu',num2str(mu),'_Nee',num2str(Nee),'_new.ds'),ds,dat_out,'fT');
    disp('Remember to unselect 3rd order gradiometer correction back in DataEditor!!')
end
if outputflag==2  %% save cleaned session out    
    meg.data=data_dssp;     % precision is now double 
    session_name(end-3:end) = [];
    save([session_name '_dssp_Mu_' num2str(mu) '_Nee_' num2str(Nee) '.mat'], 'coreg','fig',...
        'meg','Lp','voxels','voxelsize','sessionfile'); 
end   
%
figure();
%subplot(2,1,1);plot(dat_in(:,32:end,1));title('before DSSP'); %works only for ctfds
%subplot(2,1,2);plot(dat_out(:,32:end,1));title('after DSSP');
%
subplot(2,1,1);plot(dat_in_meg(:,:,1));title('before DSSP');
subplot(2,1,2);plot(data_dssp(:,:,1));title('after DSSP');
