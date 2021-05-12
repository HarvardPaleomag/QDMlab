function [nFolders] = automatic_input_ui__(nFolders, kwargs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
% Parameters
% ----------
%     nFolders: cell
%     type: str ['dir']
%     MultiSelect: str ['off']
%     single: bool [false]
%
% Returns
% -------
%   nFolders: cell

arguments
    nFolders
    kwargs.type = 'dir';
    kwargs.title = 'select file/directory';
    kwargs.MultiSelect = 'off';
    kwargs.single = false;
end

if strcmp(nFolders, 'none')
	nFolders = {};
    logMsg('input', kwargs.title, 1,0);

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

% check if any file was selected
if strcmp(nFolders{1}(2), '/')
    error('<>   ERROR: NO files/folders selected, please specify or pick files/folders.')
end

if kwargs.single
    nFolders = nFolders{:};
end

end

