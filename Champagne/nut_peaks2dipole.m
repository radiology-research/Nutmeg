function nut_peaks2dipole(peaks_file, fcm_dip_file, dip_header_file)
%function nut_peaks2dipole(peaks_file, fcm_dip_file, dip_header_file)

%peaks_file:  file with peaks.latency (nx1), 
%                   peaks.coords (nx3), peaks.amp (nx1), peaks_label,
%                   peaks_dataset
%dip_header_file:   name 'fit.dip'--this is a copy of a dipole file with
%                   headmodel and fiducials saved in it, with the 'Dipoles' line and
%                   everything after it deleted. If you name this 'fit.dip' and put it in the
%                   same directory where this program is run, this argument is optional.

if(nargin == 2)
    dip_header_file = 'fit.dip';
    if exist('fit.dip','file')~=2
        error('No .dip header file specified; fit.dip not found in current directory');
    end
end

load(peaks_file)


peaks_count = size(peaks.latency,1);

colour = 3*ones(peaks_count,1) %cyan
shape  = 0*ones(peaks_count,1)  %filled circle
if length(peaks.label) == 0
    peaks.label = 'nutmeg_peaks';
end
if length(peaks.dataset) == 0
    peaks.dataset = 'ctf.ds';
end

for jj = 1:peaks_count
    dipole_flag_index(jj) = jj;
end

dipole_flags  = [colour shape ones(peaks_count,3)];


fileID = fopen(fcm_dip_file,'a');
fprintf(fileID,'Dipoles\n');
fprintf(fileID,'{\n');
fprintf(fileID,'        // Dipole parameters ...\n');
fprintf(fileID,'        // xp (cm)         yp (cm)       zp (cm)    xo        yo         zo        Mom(nAm)        Label\n');
for jj = 1:peaks_count
    fprintf(fileID,'   \t%d: %.1f      %.1f      %.1f    0     0    0     %.2f         %s  \n',jj,peaks.coords(jj,:),peaks.amp(jj,:),peaks.label);
end
fprintf(fileID,'}\n');

fprintf(fileID,'\n');
fprintf(fileID,'Dipole_Flags\n');
fprintf(fileID,'{\n');
fprintf(fileID,'\t//Colour\tShape\tshow dir.\tshowlabel\tshow err\n');
for jj = 1:peaks_count
    fprintf(fileID,'     \t%d:    %d    %d    %d    %d    %d\n',jj,dipole_flags(jj,:));
end
fprintf(fileID,'}\n');

fprintf(fileID,'\n');
fprintf(fileID,'Dipole_FitInfo\n');
fprintf(fileID,'{\n');
fprintf(fileID,'\t//\tTrial     Start    Latency (s)    N_Pts   Error    MEG_Error     EEG_Error     ErrType  FitType  Good\tDataset\n');
for jj = 1:peaks_count
    fprintf(fileID,'     \t%d:    1         1         %.5f      %d      %d     %d     %d     %d      %d      %d     %s \n',jj,peaks.latency(jj)/1000,[1 0 0 0 0 0 1],peaks.dataset);
end
fprintf(fileID,'}\n');
fclose(fileID);

output_fitfcm_file = ['fit_' peaks.label '.dip'];  %name output file fit_cluster1p0.05.dip for example
system(['cat ' dip_header_file ' ' fcm_dip_file ' > ' output_fitfcm_file]) %concatenate fit.dip with cluster1p0.05.dip