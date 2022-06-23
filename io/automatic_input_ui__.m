function [input] = automatic_input_ui__(input, kwargs)
%[input] = automatic_input_ui__(input; 'type', 'title', 'MultiSelect', 'single', 'filter='*')
% Lets you pick a file/folder if
%
% Parameters
% ----------
%     input: cell
%     type: str ['dir']
%     MultiSelect: str ['off']
%     single: bool [false]
%
% Returns
% -------
%   input: cell

arguments
    input
    kwargs.type = 'dir';
    kwargs.title = 'select file/directory';
    kwargs.MultiSelect = 'off';
    kwargs.single = false;
    kwargs.filter='*.mat'
end

if isstruct(input)
    return
end

if strcmp(input, 'none')
	input = {};
    logMsg('input', kwargs.title, 1,0);

    if strcmp(kwargs.type, 'dir')
        
        if strcmp(kwargs.MultiSelect,'on')
            path = 1;
            while path ~= 0
                path = uigetdir('title', [kwargs.title ' cancel to stop selecting']);
                if path ~= 0
                    input{end+1} = path;
                end
            end
        else
           input = {uigetdir('title', kwargs.title)};
        end
        
    end
    if strcmp(kwargs.type, 'file')
        [file,path] = uigetfile(kwargs.filter, kwargs.title, 'MultiSelect', kwargs.MultiSelect);
        if strcmp(kwargs.MultiSelect,'on')
            for f = file
                f = fullfile(path, f);
                input{end+1} = f{:};
            end
        else
            input = {fullfile(path, file)};
        end
    end
end

input = correct_cell_shape(input);

% check if any file was selected
if size(input,1) == 1
    if strcmp(input{1}, '/')
        error('<>   ERROR: NO files/folders selected, please specify or pick files/folders.')
    end
end

if kwargs.single
    input = input{:};
end

end

