function new_DSSP_CTF(foldername, filename,subjectnumber)

Still need to add in the change on line 144 for WriteCTF

%addpath('/data/research_meg11/ccai/DSSP_Champ_ctf/ctf_fieldtrip');
%addpath('/data/research_meg11/ccai/DSSP_Champ_ctf/DSSP_new');
%addpath('/data/research_meg11/ccai/DSSP_Champ_ctf/ctf_code');
%addpath('/data/research_meg11/ccai/DSSP_Champ_ctf');
%addpath(genpath('/data/research_meg10/ccai/ccai/public_data/fieldtrip/fieldtrip-20160222'));

folder = pwd;
pattern = '*session*';
sessionfiles = dir(fullfile(folder, pattern));
pattern2='C*_*_*';
files2 = dir(fullfile(folder, pattern2));

for dsnumber = 1:length(sessionfiles)
    %        foldername = 'C3192A_EPI_20181205_09-LD2.ds'
    foldername = files2(dsnumber).name;
    filename = sessionfiles(dsnumber).name;
    
    %foldername = ['C4319A_SEF_20221212_0' num2str(dsnumber) '-RLipAv-filt1to40Hz.ds']

    which writeCTFds
    mu = 50;
    Nee = 10;
    
    
    % read .ds file
    ds = readCTFds(foldername)
    
    % get MEG channel index
    chan=cellstr(ds.res4.chanNames);
    j=0;
    for i=1:size(chan)
        if strcmp('M',chan{i}(1))
            j = j+1;
            ind(j)=i;
        end
    end
    % usually MEG channel index = 32:302
    
    %dcoffset
    %filter
    %nutmeg
    
    
    % get CTF data
    dat_in = getCTFdata(ds); %, [1;36000;1], 1:323, 'T', 'double');
    dat_in_meg = dat_in(:,ind,:);
    
    %remove dc offset
    for i=1:size(dat_in_meg,2)
        for j = 1:size(dat_in_meg,3)
            dat_in_meg(:,i,j) = detrend(dat_in_meg(:,i,j),'linear');%-mean(dat_in_meg(:,i,j));
        end
    end
    
    % load leadfiled information and MEG data
    % [filename,folder] = uigetfile('*.mat','Select MEG data file')
    load(fullfile(folder, filename));
    % load('/data/research_meg10/jjxu_t/C2621Acm10/C2621Acm10_01.mat');
    
    % prepare file to log optimal Nee for future reference
    [~,name,~]=fileparts(strcat(folder,filename));
    fid = fopen(strcat(folder,name,'_opNee.txt'),'w+');
    fprintf(fid, 'This file logs the opitimal Nee determined by its math definition in algorithm for furture reference.\n6000 data points per DSSP treatment\n');
    
    % run DSSP
    data_dssp = [];
    
    Fc = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));
    [nt,nk,ntrial]=size(dat_in_meg);
    
    
    
    for i=1:ntrial
        str = [' trial = ',num2str(i)];
        disp(str);
        [tmp,Nee_op] = MYptsss(dat_in_meg(:,:,i)',Fc, mu, Nee);
        data_dssp(:,:,i) = tmp';
        fprintf(fid,'trial %d, Nee_op = %d\n',i,Nee_op);
        close all;
    end
    Nee = Nee_op;
    %
    %
    %
    % if nt<10000
    %     for i=1:ntrial
    %         str = [' trial = ',num2str(i)];
    %         disp(str);
    %         [tmp,Nee_op] = MYptsss(dat_in_meg(:,:,i)',Fc, mu, Nee);
    %         data_dssp(:,:,i) = tmp';
    %         fprintf(fid,'trial %d, Nee_op = %d\n',i,Nee_op);
    %     end
    % else
    %     temp = round(size(dat_in_meg,1)/6000);
    %     for i=1:ntrial
    %     str = [' trial = ',num2str(i)];
    %     disp(str);
    %         for j= 1:temp-1
    %                 [tmp,Nee_op] = MYptsss(dat_in_meg((j-1)*6000+1:j*6000,:,i)',Fc, mu, Nee);
    %                  data_dssp((j-1)*6000+1:j*6000,:,i) = tmp';
    %                 fprintf(fid,'trial %d, Nee_op = %d\n',i,Nee_op);
    %         end
    %         [tmp,Nee_op] = MYptsss(dat_in_meg((temp-1)*6000+1:nt,:,i)',Fc, mu, Nee);
    %          data_dssp((temp-1)*6000+1:nt,:,i) = tmp';
    %         fprintf(fid,'trial %d, Nee_op = %d\n',i,Nee_op);
    %     end
    % end
    
    % close optimal Nee logfile
    fclose(fid);
    
    
    % rewrite DSSP-cleaned data into CTF file
    dat_out = dat_in;
    % clear dat_in;
    dat_out(:,ind,:)=data_dssp(:,:,:);
    
    %
    name = ds.baseName;
    ds_temp=writeCTFds(strcat(folder,name,'_Dssp','_mu',num2str(mu),'_Nee',num2str(Nee),'_new.ds'),ds,dat_out,'fT');
    
    % cd ctf_code/
    % writeCTFds(strcat(folder,name,'_Dssp','_mu',num2str(mu),'_Nee',num2str(Nee),'_2012.ds'),ds,dat_out,'fT');
    
    % ds_temp=writeCTFds(strcat(folder,name,'_tsss_Lin',num2str(Linmax),'_Lout',num2str(Loutmax),'.ds'),ds,dat_out,'fT');
    
    % % check to make sure the data write is correct
    % foldername2 = uigetdir();
    % ds2 = readCTFds(foldername2)
    % dat2 = getCTFdata(ds2);
    %
    % figure;
    % subplot(4,1,1); plot(dat_in(:,32:302,1));
    % subplot(4,1,2); plot(data_dssp(:,:,1));
    % subplot(4,1,3); plot(dat_out(:,32:302,1));
    % subplot(4,1,4); plot(dat2(:,32:302,1));
    
    figure();
    subplot(2,1,1);plot(dat_in(:,32:end,1));title('before DSSP');
    subplot(2,1,2);plot(dat_out(:,32:end,1));title('after DSSP');
    %
    
    curpath = pwd;
    % cd(folder);
    % save([name,'_dssp_Nee',num2str(Nee)],'dat_in','dat_out','-v7.3');
    % cd(curpath);
    
    
    % uisave({'dat_in','dat_out'},strcat(name,'_dssp_Nee',num2str(Nee)));
    
end
