function runimcoh(rawfile,chfile,ftype,lofrq,hifrq,n,stp,nfft)
%RUNIMCOH  Imaginary coherency of a 2D signal.
%
%  runimcoh(rawfile,chfile,ftype,lofrq,hifrq)

if (nargin==0 || strcmp(rawfile,'-help'))
    fprintf('\nusage: runimcoh rawfile chfile ftype lofrq hifrq\n\n')
    return
end

fprintf('\nRUNIMCOH   Calculates imaginary coherency across time.\n')
fprintf(  '------------------------------------------------------\n')
fprintf('runimcoh %s %s %s %s %s\n',rawfile,chfile,ftype,lofrq,hifrq)

if nargin<8
    nfft=256;
elseif class(nfft)==string
    nfft=str2num(nfft);
else 
    nfft=nfft
end
if (nargin<7 || isempty(stp))
    stp=1;
elseif class(stp)==string
    stp=str2num(stp);
else
    stp=stp
end
if (nargin<6 || isempty(n))
    n = 2;
elseif class(n)==string
    n=str2num(n);
else
    n=n
end 
if strcmp(class(lofrq),'double')==0
    lofrq=str2num(lofrq);
end
if strcmp(class(hifrq),'double')==0
    hifrq=str2num(hifrq);
end

load(rawfile);
clear voxels

[chpath,chfile,chext]=fileparts(chfile);
if isempty(chext), chext='.txt'; end
ch=load(fullfile(chpath,[chfile chext]));


Fs=round(1/(t(2)-t(1)));
len=size(V,1);

wdw=hanning(round(n*Fs));

n = round(n*Fs);
fsc = 1/stp;
stp = round(stp*Fs);

% Main part
frq=[max(lofrq,1/(n/Fs)) min(hifrq,Fs/2)];

nbt = floor((len-n)/stp) + 1;                           % number of time points
ffft=0:Fs/nfft:Fs/2;
select = find(ffft>=frq(1) & ffft<=frq(2));

ncomp=size(ch,1);

% freq vectors
if strcmpi(ftype,'b'), f = [ffft(select(1)) ffft(select(end))];
else,   f  = ffft(select);
end

if strcmpi(ftype,'b')
    Coh=zeros(1,ncomp);
else
    Coh=zeros(length(select),ncomp);
end

fprintf('Calculating Imaginary Coherency: ')
perc10=floor(ncomp/10);
if perc10>0, step10=[perc10:perc10:ncomp];
else, step10=ncomp;
end, clear perc10
for c=1:ncomp
    x = V(:,ch(c,1));
    y = V(:,ch(c,2));
    index = 1:n;

    Pxx = zeros(length(select),1); 
    Pyy = zeros(length(select),1); 
    Pxy = zeros(length(select),1); 

    for k=1:nbt
        xw = x(index,:);                      % Windowing
        yw = y(index,:);
        xw = wdw.*detrend(xw);
        yw = wdw.*detrend(yw);
        Xx = fft(xw,nfft);
        Yy = fft(yw,nfft);
        Xx2 = abs(Xx).^2;
        Yy2 = abs(Yy).^2;
        Xy2 = Yy.*conj(Xx);
        Pxx = Pxx + Xx2(select);
        Pyy = Pyy + Yy2(select);
        Pxy = Pxy + Xy2(select);
        index = index + stp;    
    end
    switch lower(ftype)
    case 'b'
        Pxx = sum(Pxx);
        Pyy = sum(Pyy);
        Pxy = sum(Pxy); 
        Coh(c) = imag(Pxy/sqrt(Pxx*Pyy));             % coherence function estimate
    case 'a'
        Coh(:,c) = imag(Pxy./sqrt(Pxx.*Pyy));        % coherence function estimate 
    end
    if any(step10==c)
        fprintf('...%d%%',ceil((c/ncomp)*100))
    end
end
fprintf('\n')

IC.coh=Coh;
IC.frq=f;
IC.comps=ch;
IC.time=[t(1) t(end)];

if ~exist('imcoh','dir'), mkdir imcoh, end
save(fullfile('imcoh',['I' chfile(2:end)]),'IC')

fprintf('Done (runimcoh).\n')