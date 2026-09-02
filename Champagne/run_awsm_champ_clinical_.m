folder = pwd;
pattern = '*session*';
% Get a list of files that match the search pattern
sessionfiles = dir(fullfile(folder, pattern));
% Display the names of the files

choice = questdlg('Choose an option', ...
    'Champagne Type', ...
    'Clinical','Research','Clinical');

switch choice
    case 'Clinical'
        words = ["D2.mat", "lip.mat", "AEF.mat"];
        for i=1:length(sessionfiles)
            sf=fullfile(sessionfiles(i).folder,sessionfiles(i).name);
            for j=1:length(words)
                if contains(sessionfiles(i).name,words(j),IgnoreCase=true)
                    switch words(j)
                        case 'D2.mat'
                            champagne_wrapper(sf,-125,-5,0,150,0);
                        case 'lip.mat'
                            champagne_wrapper(sf,-125,-5,0,150,0);
                        case 'AEF.mat'
                            champagne_wrapper(sf,-100,-5,0,200,0);
                    end
                end
            end
        end
    case 'Research'
        disp('You chose Option 2.')
        [filename, filepath] = uigetfile('session*', 'Select a session file...');
        if isequal(filename, 0) || isequal(filepath, 0)
            disp('No file selected.');
        else
            disp(['Selected file: ' fullfile(filepath, filename)]);
        end
        sf=fullfile(filepath, filename);
        prompt = {'Starting Control Window:', 'Ending Control Windo:', ...
            'Starting Active Window:','Ending Active Window:','Optional VCS (Keep as 0 outside of tuning)'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'-125','-5','0','150','0'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        champagne_wrapper(sf,str2num(answer{1}),str2num(answer{2}),...
            str2num(answer{3}),str2num(answer{4}),str2num(answer{5}));
end
