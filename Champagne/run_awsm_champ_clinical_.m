% Uncomment the line with desired response localization type and update path
% Typically, will not have to change values
% WARNING: Make sure to change the first value (control window) for SEF to -100 for Epilepsy cases!!!


%Change this so it doesnt call nutmeg again

folder = pwd;
pattern = '*session*';
% Get a list of files that match the search pattern
sessionfiles = dir(fullfile(folder, pattern));
% Display the names of the files
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