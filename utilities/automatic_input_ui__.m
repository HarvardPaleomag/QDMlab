function [nFolders] = automatic_input_ui__(nFolders, kwargs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    nFolders
    kwargs.type = 'dir';
    kwargs.title = 'select file/directory';
    kwargs.MultiSelect = 'off';
end

if nFolders == 'none'
	nFolders = {};

    if strcmp(kwargs.type, 'dir')
        
        if strcmp(kwargs.MultiSelect,'on')
            path = 1;
            while path ~= 0
                path = uigetdir('title', [kwargs.title ' cancel to stop selecting']);
                if path ~= 0
                    nFolders{end+1} = path;
                end
            end
        else
           nFolders = {uigetdir('title', kwargs.title)};
        end
        
    end
    if strcmp(kwargs.type, 'file')
        [file,path] = uigetfile('title', kwargs.title, 'MultiSelect', kwargs.MultiSelect);
        if strcmp(kwargs.MultiSelect,'on')
            for f = file
                f = fullfile(path, f);
                nFolders{end+1} = f{:};
            end
        else
            nFolders = {fullfile(path, file)};
        end
    end
end

nFolders = correct_cell_shape(nFolders);
end
